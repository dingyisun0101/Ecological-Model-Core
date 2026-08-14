//! Ecological interaction matrices, reproducible recipes, and verified artifacts.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use physics_in_parallel::math::prelude::{DenseMatrix, MatrixError};
use physics_in_parallel::rng::{IndexedRng, RngConfig, RngConfigError};
use scientific_workflow::artifact::{
    ArtifactDescriptor, ArtifactDisposition, ArtifactError, ArtifactLoadError,
    load_verified_artifact, persist_artifact,
};
use scientific_workflow::execution::ExecutionScope;
use scientific_workflow::rng_record::{RngRecord, RngRecordError};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use thiserror::Error;

pub const INTERACTION_MATRIX_FORMAT: &str = "ecological.interaction-matrix.v1";
pub const INTERACTION_MATRIX_METADATA_KEY: &str = "interaction_matrix";
pub const INTERACTION_GENERATOR_RNG_NAMESPACE: &str = "ecological_model_core.interaction_matrix";
pub const INTERACTION_GENERATOR_IDENTITY: &str = "ecological_model_core.interaction_matrix";
pub const INTERACTION_GENERATOR_VERSION: &str = "1";

const DOMAIN_CONNECTANCE: u64 = 0x5c1a_9f20_f678_314d;
const DOMAIN_FIRST_NORMAL: u64 = 0x9841_d60a_334b_c8e7;
const DOMAIN_SECOND_NORMAL: u64 = 0xa72e_1b49_963c_05fd;

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum MatrixNormalization {
    None,
    SqrtSpecies,
}

impl MatrixNormalization {
    fn divisor(self, species: usize) -> f64 {
        match self {
            Self::None => 1.0,
            Self::SqrtSpecies => (species as f64).sqrt(),
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "kind", content = "value", rename_all = "snake_case")]
pub enum DiagonalPolicy {
    Zero,
    Constant(f64),
}

impl DiagonalPolicy {
    fn value(self) -> f64 {
        match self {
            Self::Zero => 0.0,
            Self::Constant(value) => value,
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SignStructure {
    Competition,
    Mutualism,
    ConsumerResource,
}

/// Reproducible scientific instructions for an ecological matrix realization.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "family", rename_all = "snake_case", deny_unknown_fields)]
pub enum InteractionMatrixRecipe {
    /// Exact Dispatcher family: zero diagonal and `A[j,i] = -A[i,j]`.
    AntisymmetricGaussian {
        scale: f64,
        #[serde(default = "sqrt_species")]
        normalization: MatrixNormalization,
        #[serde(default)]
        rng: RngConfig,
    },
    /// Independent directed Gaussian coefficients with optional sparsity.
    IndependentGaussian {
        mean: f64,
        standard_deviation: f64,
        #[serde(default = "one")]
        connectance: f64,
        #[serde(default = "zero_diagonal")]
        diagonal: DiagonalPolicy,
        #[serde(default = "sqrt_species")]
        normalization: MatrixNormalization,
        #[serde(default)]
        rng: RngConfig,
    },
    /// Gaussian reciprocal pairs with Pearson correlation in `[-1,1]`.
    CorrelatedGaussian {
        mean: f64,
        standard_deviation: f64,
        reciprocal_correlation: f64,
        #[serde(default = "one")]
        connectance: f64,
        #[serde(default = "zero_diagonal")]
        diagonal: DiagonalPolicy,
        #[serde(default = "sqrt_species")]
        normalization: MatrixNormalization,
        #[serde(default)]
        rng: RngConfig,
    },
    /// Pairwise signs constrained to a common ecological relationship class.
    SignStructuredGaussian {
        structure: SignStructure,
        scale: f64,
        #[serde(default = "one")]
        connectance: f64,
        #[serde(default = "zero_diagonal")]
        diagonal: DiagonalPolicy,
        #[serde(default = "sqrt_species")]
        normalization: MatrixNormalization,
        #[serde(default)]
        rng: RngConfig,
    },
}

const fn sqrt_species() -> MatrixNormalization {
    MatrixNormalization::SqrtSpecies
}
const fn zero_diagonal() -> DiagonalPolicy {
    DiagonalPolicy::Zero
}
const fn one() -> f64 {
    1.0
}

impl InteractionMatrixRecipe {
    pub fn validate(&self, species: usize) -> Result<(), InteractionRecipeError> {
        if species == 0 {
            return Err(InteractionRecipeError::EmptySpecies);
        }
        match self {
            Self::AntisymmetricGaussian { scale, .. } => {
                require_nonnegative_finite("scale", *scale)?;
            }
            Self::IndependentGaussian {
                mean,
                standard_deviation,
                connectance,
                diagonal,
                ..
            }
            | Self::CorrelatedGaussian {
                mean,
                standard_deviation,
                connectance,
                diagonal,
                ..
            } => {
                require_finite("mean", *mean)?;
                require_nonnegative_finite("standard_deviation", *standard_deviation)?;
                require_probability(*connectance)?;
                require_finite("diagonal", diagonal.value())?;
            }
            Self::SignStructuredGaussian {
                scale,
                connectance,
                diagonal,
                ..
            } => {
                require_nonnegative_finite("scale", *scale)?;
                require_probability(*connectance)?;
                require_finite("diagonal", diagonal.value())?;
            }
        }
        if let Self::CorrelatedGaussian {
            reciprocal_correlation,
            ..
        } = self
            && (!reciprocal_correlation.is_finite()
                || !(-1.0..=1.0).contains(reciprocal_correlation))
        {
            return Err(InteractionRecipeError::InvalidParameter {
                name: "reciprocal_correlation",
                value: *reciprocal_correlation,
            });
        }
        Ok(())
    }

    pub fn generate(&self, species: usize) -> Result<InteractionMatrix, InteractionRecipeError> {
        self.validate(species)?;
        let rng = IndexedRng::new(self.rng())?;
        let resolved_recipe = self.with_rng(rng.rng_config());
        let mut values = vec![0.0; species * species];
        match &resolved_recipe {
            Self::AntisymmetricGaussian {
                scale,
                normalization,
                ..
            } => {
                let divisor = normalization.divisor(species);
                for row in 0..species {
                    for column in (row + 1)..species {
                        // Preserve Dispatcher's established realization exactly: its indexed
                        // coordinates are `(0, species, row, column)`.
                        let value =
                            rng.standard_normal(0, species as u64, row as u64, column as u64)
                                * scale
                                / divisor;
                        values[row * species + column] = value;
                        values[column * species + row] = -value;
                    }
                }
            }
            Self::IndependentGaussian {
                mean,
                standard_deviation,
                connectance,
                diagonal,
                normalization,
                ..
            } => {
                let divisor = normalization.divisor(species);
                for row in 0..species {
                    for column in 0..species {
                        values[row * species + column] = if row == column {
                            diagonal.value()
                        } else if connected(rng, species, row, column, *connectance, false) {
                            (mean
                                + standard_deviation
                                    * normal(rng, DOMAIN_FIRST_NORMAL, species, row, column, 0))
                                / divisor
                        } else {
                            0.0
                        };
                    }
                }
            }
            Self::CorrelatedGaussian {
                mean,
                standard_deviation,
                reciprocal_correlation,
                connectance,
                diagonal,
                normalization,
                ..
            } => {
                fill_diagonal(&mut values, species, diagonal.value());
                let divisor = normalization.divisor(species);
                let independent_weight = (1.0 - reciprocal_correlation.powi(2)).sqrt();
                for row in 0..species {
                    for column in (row + 1)..species {
                        if !connected(rng, species, row, column, *connectance, true) {
                            continue;
                        }
                        let first = normal(rng, DOMAIN_FIRST_NORMAL, species, row, column, 0);
                        let second = normal(rng, DOMAIN_SECOND_NORMAL, species, row, column, 1);
                        values[row * species + column] =
                            (mean + standard_deviation * first) / divisor;
                        values[column * species + row] = (mean
                            + standard_deviation
                                * (reciprocal_correlation * first + independent_weight * second))
                            / divisor;
                    }
                }
            }
            Self::SignStructuredGaussian {
                structure,
                scale,
                connectance,
                diagonal,
                normalization,
                ..
            } => {
                fill_diagonal(&mut values, species, diagonal.value());
                let magnitude = scale / normalization.divisor(species);
                for row in 0..species {
                    for column in (row + 1)..species {
                        if !connected(rng, species, row, column, *connectance, true) {
                            continue;
                        }
                        let first = normal(rng, DOMAIN_FIRST_NORMAL, species, row, column, 0).abs()
                            * magnitude;
                        let second = normal(rng, DOMAIN_SECOND_NORMAL, species, row, column, 1)
                            .abs()
                            * magnitude;
                        let (forward, reverse) = match structure {
                            SignStructure::Competition => (-first, -second),
                            SignStructure::Mutualism => (first, second),
                            SignStructure::ConsumerResource => (first, -second),
                        };
                        values[row * species + column] = forward;
                        values[column * species + row] = reverse;
                    }
                }
            }
        }
        let generator = GeneratorProvenance::new(
            INTERACTION_GENERATOR_IDENTITY,
            INTERACTION_GENERATOR_VERSION,
            serde_json::to_value(&resolved_recipe)?,
            Some(rng.rng_config()),
        )?;
        Ok(InteractionMatrix::from_generated(
            DenseMatrix::try_from_vec(species, species, values)?,
            species,
            generator,
        )?)
    }

    pub const fn rng(&self) -> RngConfig {
        match self {
            Self::AntisymmetricGaussian { rng, .. }
            | Self::IndependentGaussian { rng, .. }
            | Self::CorrelatedGaussian { rng, .. }
            | Self::SignStructuredGaussian { rng, .. } => *rng,
        }
    }

    fn with_rng(&self, resolved: RngConfig) -> Self {
        let mut recipe = self.clone();
        match &mut recipe {
            Self::AntisymmetricGaussian { rng, .. }
            | Self::IndependentGaussian { rng, .. }
            | Self::CorrelatedGaussian { rng, .. }
            | Self::SignStructuredGaussian { rng, .. } => *rng = resolved,
        }
        recipe
    }
}

fn fill_diagonal(values: &mut [f64], species: usize, diagonal: f64) {
    for index in 0..species {
        values[index * species + index] = diagonal;
    }
}

fn normal(
    rng: IndexedRng,
    domain: u64,
    species: usize,
    row: usize,
    column: usize,
    component: u64,
) -> f64 {
    rng.standard_normal(
        0,
        domain ^ species as u64,
        row as u64,
        (column as u64) ^ component,
    )
}

fn connected(
    rng: IndexedRng,
    species: usize,
    row: usize,
    column: usize,
    connectance: f64,
    unordered: bool,
) -> bool {
    if connectance >= 1.0 {
        return true;
    }
    let (row, column) = if unordered && row > column {
        (column, row)
    } else {
        (row, column)
    };
    rng.unit_f64(
        0,
        DOMAIN_CONNECTANCE ^ species as u64,
        row as u64,
        column as u64,
        0,
    ) < connectance
}

fn require_finite(name: &'static str, value: f64) -> Result<(), InteractionRecipeError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(InteractionRecipeError::InvalidParameter { name, value })
    }
}

fn require_nonnegative_finite(
    name: &'static str,
    value: f64,
) -> Result<(), InteractionRecipeError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(InteractionRecipeError::InvalidParameter { name, value })
    }
}

fn require_probability(value: f64) -> Result<(), InteractionRecipeError> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(())
    } else {
        Err(InteractionRecipeError::InvalidParameter {
            name: "connectance",
            value,
        })
    }
}

#[derive(Clone, Debug)]
pub struct InteractionMatrix {
    values: Arc<DenseMatrix<f64>>,
    provenance: InteractionProvenance,
}

impl InteractionMatrix {
    pub fn from_matrix(
        values: DenseMatrix<f64>,
        species: usize,
    ) -> Result<Self, InteractionMatrixError> {
        Self::resolve(
            Arc::new(values),
            species,
            InteractionProvenance::InMemory { label: None },
        )
    }

    pub fn from_shared(
        values: Arc<DenseMatrix<f64>>,
        species: usize,
    ) -> Result<Self, InteractionMatrixError> {
        Self::resolve(
            values,
            species,
            InteractionProvenance::InMemory { label: None },
        )
    }

    pub fn from_labeled_matrix(
        values: DenseMatrix<f64>,
        species: usize,
        label: impl Into<String>,
    ) -> Result<Self, InteractionMatrixError> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(InteractionMatrixError::EmptyLabel);
        }
        Self::resolve(
            Arc::new(values),
            species,
            InteractionProvenance::InMemory { label: Some(label) },
        )
    }

    pub fn from_rows(rows: Vec<Vec<f64>>, species: usize) -> Result<Self, InteractionMatrixError> {
        let row_count = rows.len();
        let column_count = rows.first().map_or(0, Vec::len);
        for (row, values) in rows.iter().enumerate() {
            if values.len() != column_count {
                return Err(InteractionMatrixError::RaggedRows {
                    row,
                    expected: column_count,
                    actual: values.len(),
                });
            }
        }
        let values = DenseMatrix::try_from_vec(
            row_count,
            column_count,
            rows.into_iter().flatten().collect(),
        )?;
        Self::resolve(Arc::new(values), species, InteractionProvenance::Inline)
    }

    pub fn load_json(
        path: impl Into<PathBuf>,
        species: usize,
    ) -> Result<Self, InteractionMatrixError> {
        let path = path.into();
        let bytes = fs::read(&path).map_err(|source| InteractionMatrixError::Io {
            path: path.clone(),
            source,
        })?;
        Self::from_json_bytes(bytes, path, species, None)
    }

    pub fn from_generated(
        values: DenseMatrix<f64>,
        species: usize,
        generator: GeneratorProvenance,
    ) -> Result<Self, InteractionMatrixError> {
        Self::resolve(
            Arc::new(values),
            species,
            InteractionProvenance::Generated { generator },
        )
    }

    pub fn generate(
        species: usize,
        recipe: &InteractionMatrixRecipe,
    ) -> Result<Self, InteractionRecipeError> {
        recipe.generate(species)
    }

    pub fn species(&self) -> usize {
        self.values.rows()
    }
    pub fn values(&self) -> &DenseMatrix<f64> {
        &self.values
    }
    pub fn shared_values(&self) -> Arc<DenseMatrix<f64>> {
        Arc::clone(&self.values)
    }
    #[inline]
    pub fn coefficient(&self, row: usize, column: usize) -> f64 {
        self.values.get(row as isize, column as isize)
    }
    #[inline]
    pub fn mul_vector_into(&self, input: &[f64], output: &mut [f64]) -> Result<(), MatrixError> {
        self.values.mul_vector_into(input, output)
    }
    pub const fn provenance(&self) -> &InteractionProvenance {
        &self.provenance
    }
    pub fn generator_rng_record(&self) -> Result<Option<RngRecord>, RngRecordError> {
        self.provenance.generator_rng_record()
    }

    fn from_json_bytes(
        bytes: Vec<u8>,
        path: PathBuf,
        species: usize,
        generator: Option<GeneratorProvenance>,
    ) -> Result<Self, InteractionMatrixError> {
        let values =
            serde_json::from_slice(&bytes).map_err(|source| InteractionMatrixError::Json {
                path: path.clone(),
                source,
            })?;
        let provenance = generator.map_or(InteractionProvenance::JsonFile { path }, |generator| {
            InteractionProvenance::Generated { generator }
        });
        Self::resolve(Arc::new(values), species, provenance)
    }

    fn resolve(
        values: Arc<DenseMatrix<f64>>,
        species: usize,
        provenance: InteractionProvenance,
    ) -> Result<Self, InteractionMatrixError> {
        if species == 0 {
            return Err(InteractionMatrixError::EmptySpecies);
        }
        let rows = values.rows();
        let columns = values.cols();
        if rows != columns {
            return Err(InteractionMatrixError::NonSquare { rows, columns });
        }
        if rows != species {
            return Err(InteractionMatrixError::SpeciesMismatch {
                expected: species,
                actual: rows,
            });
        }
        for row in 0..rows {
            for column in 0..columns {
                let value = values.get(row as isize, column as isize);
                if !value.is_finite() {
                    return Err(InteractionMatrixError::NonFiniteEntry { row, column, value });
                }
            }
        }
        Ok(Self { values, provenance })
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum InteractionSourceKind {
    InMemory,
    Inline,
    JsonFile,
    Generated,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum InteractionProvenance {
    InMemory { label: Option<String> },
    Inline,
    JsonFile { path: PathBuf },
    Generated { generator: GeneratorProvenance },
}

impl InteractionProvenance {
    pub const fn kind(&self) -> InteractionSourceKind {
        match self {
            Self::InMemory { .. } => InteractionSourceKind::InMemory,
            Self::Inline => InteractionSourceKind::Inline,
            Self::JsonFile { .. } => InteractionSourceKind::JsonFile,
            Self::Generated { .. } => InteractionSourceKind::Generated,
        }
    }
    pub const fn generator(&self) -> Option<&GeneratorProvenance> {
        match self {
            Self::Generated { generator } => Some(generator),
            _ => None,
        }
    }
    pub fn generator_rng_record(&self) -> Result<Option<RngRecord>, RngRecordError> {
        self.generator()
            .map_or(Ok(None), GeneratorProvenance::rng_record)
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GeneratorProvenance {
    identity: String,
    version: String,
    recipe: Value,
    rng: Option<RngConfig>,
}

impl GeneratorProvenance {
    pub fn new(
        identity: impl Into<String>,
        version: impl Into<String>,
        recipe: Value,
        rng: Option<RngConfig>,
    ) -> Result<Self, InteractionMatrixError> {
        let identity = identity.into();
        let version = version.into();
        if identity.trim().is_empty() {
            return Err(InteractionMatrixError::InvalidGeneratorLabel { field: "identity" });
        }
        if version.trim().is_empty() {
            return Err(InteractionMatrixError::InvalidGeneratorLabel { field: "version" });
        }
        if rng.is_some_and(|value| value.seed().is_none() || value.method().is_none()) {
            return Err(InteractionMatrixError::UnresolvedGeneratorRng { identity });
        }
        Ok(Self {
            identity,
            version,
            recipe,
            rng,
        })
    }
    pub fn identity(&self) -> &str {
        &self.identity
    }
    pub fn version(&self) -> &str {
        &self.version
    }
    pub const fn recipe(&self) -> &Value {
        &self.recipe
    }
    pub const fn rng(&self) -> Option<RngConfig> {
        self.rng
    }
    pub fn rng_record(&self) -> Result<Option<RngRecord>, RngRecordError> {
        let Some(rng) = self.rng else {
            return Ok(None);
        };
        let method = rng.method().expect("generator RNG is resolved");
        let mut parameters = Map::new();
        parameters.insert("recipe".to_owned(), self.recipe.clone());
        if let Some(streams) = rng.parallel_streams() {
            parameters.insert("parallel_streams".to_owned(), Value::from(streams.get()));
        }
        Ok(Some(RngRecord::new(
            INTERACTION_GENERATOR_RNG_NAMESPACE,
            format!("{}+{}", self.identity, method.name()),
            format!("{}+{}", self.version, method.version()),
            method.seed_encoding(),
            rng.encode_seed().expect("generator RNG seed is resolved"),
            Some(parameters),
        )?))
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct InteractionArtifactDescriptor {
    format: String,
    species: usize,
    #[serde(flatten)]
    artifact: ArtifactDescriptor,
    source_kind: InteractionSourceKind,
    #[serde(skip_serializing_if = "Option::is_none")]
    generator: Option<GeneratorProvenance>,
}

impl InteractionArtifactDescriptor {
    pub fn format(&self) -> &str {
        &self.format
    }
    pub const fn species(&self) -> usize {
        self.species
    }
    pub const fn shape(&self) -> [usize; 2] {
        [self.species, self.species]
    }
    pub fn sha256(&self) -> &str {
        self.artifact.sha256()
    }
    pub fn path(&self) -> &str {
        self.artifact.path()
    }
    pub const fn source_kind(&self) -> InteractionSourceKind {
        self.source_kind
    }
    pub const fn generator(&self) -> Option<&GeneratorProvenance> {
        self.generator.as_ref()
    }
    pub fn insert_into_metadata(&self, metadata: &mut Map<String, Value>) -> Option<Value> {
        metadata.insert(
            INTERACTION_MATRIX_METADATA_KEY.to_owned(),
            serde_json::to_value(self).expect("interaction descriptor is JSON-compatible"),
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PersistedInteraction {
    descriptor: InteractionArtifactDescriptor,
    disposition: ArtifactDisposition,
}

impl PersistedInteraction {
    pub const fn descriptor(&self) -> &InteractionArtifactDescriptor {
        &self.descriptor
    }
    pub const fn disposition(&self) -> ArtifactDisposition {
        self.disposition
    }
    pub fn into_descriptor(self) -> InteractionArtifactDescriptor {
        self.descriptor
    }
}

pub fn persist_interaction_matrix(
    scope: &ExecutionScope,
    matrix: &InteractionMatrix,
) -> Result<PersistedInteraction, InteractionArtifactError> {
    let bytes = serde_json::to_vec(matrix.values())?;
    let persisted = persist_artifact(scope, "interaction", "json", &bytes)?;
    Ok(PersistedInteraction {
        descriptor: InteractionArtifactDescriptor {
            format: INTERACTION_MATRIX_FORMAT.to_owned(),
            species: matrix.species(),
            artifact: persisted.descriptor().clone(),
            source_kind: matrix.provenance().kind(),
            generator: matrix.provenance().generator().cloned(),
        },
        disposition: persisted.disposition(),
    })
}

pub fn load_verified_interaction_matrix(
    execution_directory: impl AsRef<Path>,
    descriptor: &InteractionArtifactDescriptor,
) -> Result<InteractionMatrix, InteractionArtifactLoadError> {
    if descriptor.format != INTERACTION_MATRIX_FORMAT || descriptor.species == 0 {
        return Err(InteractionArtifactLoadError::InvalidDescriptor);
    }
    let verified = load_verified_artifact(execution_directory, &descriptor.artifact)?;
    let path = verified.path().to_path_buf();
    Ok(InteractionMatrix::from_json_bytes(
        verified.into_bytes(),
        path,
        descriptor.species,
        descriptor.generator.clone(),
    )?)
}

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum InteractionRecipeError {
    #[error("interaction recipe requires at least one species")]
    EmptySpecies,
    #[error("interaction recipe parameter {name} is invalid: {value}")]
    InvalidParameter { name: &'static str, value: f64 },
    #[error(transparent)]
    Rng(#[from] RngConfigError),
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    #[error(transparent)]
    Interaction(#[from] InteractionMatrixError),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
}

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum InteractionMatrixError {
    #[error(transparent)]
    Matrix(#[from] MatrixError),
    #[error("interaction matrix species dimension must be positive")]
    EmptySpecies,
    #[error("interaction matrix must be square, found {rows}x{columns}")]
    NonSquare { rows: usize, columns: usize },
    #[error("interaction matrix has {actual} species, expected {expected}")]
    SpeciesMismatch { expected: usize, actual: usize },
    #[error("interaction matrix row {row} has {actual} columns, expected {expected}")]
    RaggedRows {
        row: usize,
        expected: usize,
        actual: usize,
    },
    #[error("interaction matrix entry ({row}, {column}) is not finite: {value}")]
    NonFiniteEntry {
        row: usize,
        column: usize,
        value: f64,
    },
    #[error("interaction matrix label must not be empty")]
    EmptyLabel,
    #[error("failed to read interaction matrix at `{path}`")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("invalid interaction matrix JSON at `{path}`")]
    Json {
        path: PathBuf,
        #[source]
        source: serde_json::Error,
    },
    #[error("interaction generator {field} must not be empty")]
    InvalidGeneratorLabel { field: &'static str },
    #[error("interaction generator `{identity}` has unresolved RNG")]
    UnresolvedGeneratorRng { identity: String },
}

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum InteractionArtifactError {
    #[error(transparent)]
    Serialize(#[from] serde_json::Error),
    #[error(transparent)]
    Workflow(#[from] ArtifactError),
}

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum InteractionArtifactLoadError {
    #[error("invalid interaction artifact descriptor")]
    InvalidDescriptor,
    #[error(transparent)]
    Workflow(#[from] ArtifactLoadError),
    #[error(transparent)]
    Matrix(#[from] InteractionMatrixError),
}
