//! Ecological interaction matrices, reproducible recipes, and verified artifacts.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use physics_in_parallel::prelude::basic::{
    DenseMatrix, MatrixError, RandType, RngConfig, RngConfigError, TensorRandError,
    TensorRandFiller,
};
use scientific_workflow::prelude::basics::{
    ArtifactDescriptor, ArtifactDisposition, ArtifactError, ArtifactLoadError, ExecutionScope,
    RngRecord, RngRecordError, load_verified_artifact, persist_artifact,
};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use thiserror::Error;

pub const INTERACTION_ARTIFACT_REFERENCE_FORMAT: &str =
    "ecological.interaction-artifact-reference.v1";

pub const INTERACTION_MATRIX_FORMAT: &str = "ecological.interaction-matrix.v2";
pub const INTERACTION_MATRIX_METADATA_KEY: &str = "interaction_matrix";
pub const INTERACTION_GENERATOR_RNG_NAMESPACE: &str = "ecological_model_core.interaction_matrix";
pub const INTERACTION_GENERATOR_IDENTITY: &str = "ecological_model_core.interaction_matrix";
pub const INTERACTION_GENERATOR_VERSION: &str = "3";

const DOMAIN_INDEPENDENT: u64 = 0x73d8_ba6e_209f_54c1;
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
    /// Sample diagonal entries from the Gaussian distribution of the recipe.
    Sampled,
}

impl DiagonalPolicy {
    fn fixed_value(self) -> Option<f64> {
        match self {
            Self::Zero => Some(0.0),
            Self::Constant(value) => Some(value),
            Self::Sampled => None,
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
    /// Independently sample every matrix entry from a uniform distribution.
    RandomUniform {
        minimum: f64,
        maximum: f64,
        #[serde(default)]
        rng: RngConfig,
    },
    /// Independently sample every matrix entry from a Gaussian distribution.
    RandomGaussian {
        mean: f64,
        standard_deviation: f64,
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
            Self::RandomUniform {
                minimum, maximum, ..
            } => {
                require_finite("minimum", *minimum)?;
                require_finite("maximum", *maximum)?;
                if minimum >= maximum {
                    return Err(InteractionRecipeError::InvalidRange {
                        minimum: *minimum,
                        maximum: *maximum,
                    });
                }
            }
            Self::RandomGaussian {
                mean,
                standard_deviation,
                ..
            } => {
                require_finite("mean", *mean)?;
                require_nonnegative_finite("standard_deviation", *standard_deviation)?;
            }
            Self::CorrelatedGaussian {
                mean,
                standard_deviation,
                connectance,
                diagonal,
                ..
            } => {
                require_finite("mean", *mean)?;
                require_nonnegative_finite("standard_deviation", *standard_deviation)?;
                require_probability(*connectance)?;
                if let Some(value) = diagonal.fixed_value() {
                    require_finite("diagonal", value)?;
                }
            }
            Self::SignStructuredGaussian {
                scale,
                connectance,
                diagonal,
                ..
            } => {
                require_nonnegative_finite("scale", *scale)?;
                require_probability(*connectance)?;
                let value = diagonal.fixed_value().ok_or(
                    InteractionRecipeError::SampledDiagonalUnsupported {
                        family: "sign_structured_gaussian",
                    },
                )?;
                require_finite("diagonal", value)?;
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
        let kind = match self {
            Self::RandomUniform {
                minimum, maximum, ..
            } => RandType::Uniform {
                low: *minimum,
                high: *maximum,
            },
            Self::RandomGaussian {
                mean,
                standard_deviation,
                ..
            } => RandType::Normal {
                mean: *mean,
                std: *standard_deviation,
            },
            Self::CorrelatedGaussian { .. } | Self::SignStructuredGaussian { .. } => {
                RandType::Normal {
                    mean: 0.0,
                    std: 1.0,
                }
            }
        };
        let mut filler = TensorRandFiller::try_new_indexed(kind, self.rng())?;
        let resolved_recipe = self.with_rng(filler.rng_config());
        let matrix_len = species
            .checked_mul(species)
            .ok_or(MatrixError::ShapeProductOverflow {
                rows: species,
                cols: species,
            })?;
        let mut values = vec![0.0; matrix_len];
        match &resolved_recipe {
            Self::RandomUniform { .. } | Self::RandomGaussian { .. } => {
                filler.try_fill_slice_at_layout(
                    &mut values,
                    species,
                    0,
                    DOMAIN_INDEPENDENT ^ species as u64,
                )?;
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
                let samples =
                    sample_structured_inputs(&mut filler, species, matrix_len, &mut values)?;
                let first_normal = samples.first_normal;
                let second_normal = samples.second_normal;
                let divisor = normalization.divisor(species);
                for index in 0..species {
                    values[index * species + index] = diagonal.fixed_value().unwrap_or_else(|| {
                        (mean + standard_deviation * first_normal[index * species + index])
                            / divisor
                    });
                }
                let independent_weight = (1.0 - reciprocal_correlation.powi(2)).sqrt();
                for row in 0..species {
                    for column in (row + 1)..species {
                        let index = row * species + column;
                        if *connectance < 1.0 && values[index] >= *connectance {
                            continue;
                        }
                        let first = first_normal[index];
                        let second = second_normal[index];
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
                let samples =
                    sample_structured_inputs(&mut filler, species, matrix_len, &mut values)?;
                let first_normal = samples.first_normal;
                let second_normal = samples.second_normal;
                fill_diagonal(
                    &mut values,
                    species,
                    diagonal
                        .fixed_value()
                        .expect("sampled diagonal rejected during validation"),
                );
                let magnitude = scale / normalization.divisor(species);
                for row in 0..species {
                    for column in (row + 1)..species {
                        let index = row * species + column;
                        if *connectance < 1.0 && values[index] >= *connectance {
                            continue;
                        }
                        let first = first_normal[index].abs() * magnitude;
                        let second = second_normal[index].abs() * magnitude;
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
            Some(filler.rng_config()),
        )?;
        Ok(InteractionMatrix::from_generated(
            DenseMatrix::try_from_vec(species, species, values)?,
            generator,
        )?)
    }

    pub const fn rng(&self) -> RngConfig {
        match self {
            Self::RandomUniform { rng, .. }
            | Self::RandomGaussian { rng, .. }
            | Self::CorrelatedGaussian { rng, .. }
            | Self::SignStructuredGaussian { rng, .. } => *rng,
        }
    }

    fn with_rng(&self, resolved: RngConfig) -> Self {
        let mut recipe = self.clone();
        match &mut recipe {
            Self::RandomUniform { rng, .. }
            | Self::RandomGaussian { rng, .. }
            | Self::CorrelatedGaussian { rng, .. }
            | Self::SignStructuredGaussian { rng, .. } => *rng = resolved,
        }
        recipe
    }
}

struct StructuredSamples {
    first_normal: Vec<f64>,
    second_normal: Vec<f64>,
}

fn sample_structured_inputs(
    filler: &mut TensorRandFiller,
    species: usize,
    matrix_len: usize,
    output: &mut [f64],
) -> Result<StructuredSamples, TensorRandError> {
    let mut first_normal = vec![0.0; matrix_len];
    let mut second_normal = vec![0.0; matrix_len];
    filler.try_fill_slice_at_layout(
        &mut first_normal,
        species,
        0,
        DOMAIN_FIRST_NORMAL ^ species as u64,
    )?;
    filler.try_fill_slice_at_layout(
        &mut second_normal,
        species,
        0,
        DOMAIN_SECOND_NORMAL ^ species as u64,
    )?;
    filler.set_kind(RandType::Uniform {
        low: 0.0,
        high: 1.0,
    });
    filler.try_fill_slice_at_layout(output, species, 0, DOMAIN_CONNECTANCE ^ species as u64)?;
    Ok(StructuredSamples {
        first_normal,
        second_normal,
    })
}

fn fill_diagonal(values: &mut [f64], species: usize, diagonal: f64) {
    for index in 0..species {
        values[index * species + index] = diagonal;
    }
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
    pub fn from_matrix(values: DenseMatrix<f64>) -> Result<Self, InteractionMatrixError> {
        Self::resolve(
            Arc::new(values),
            InteractionProvenance::InMemory { label: None },
        )
    }

    pub fn from_shared(values: Arc<DenseMatrix<f64>>) -> Result<Self, InteractionMatrixError> {
        Self::resolve(values, InteractionProvenance::InMemory { label: None })
    }

    pub fn from_labeled_matrix(
        values: DenseMatrix<f64>,
        label: impl Into<String>,
    ) -> Result<Self, InteractionMatrixError> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(InteractionMatrixError::EmptyLabel);
        }
        Self::resolve(
            Arc::new(values),
            InteractionProvenance::InMemory { label: Some(label) },
        )
    }

    pub fn from_rows(rows: Vec<Vec<f64>>) -> Result<Self, InteractionMatrixError> {
        let row_count = rows.len();
        if row_count == 0 {
            return Err(InteractionMatrixError::EmptySpecies);
        }
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
        Self::resolve(Arc::new(values), InteractionProvenance::Inline)
    }

    pub fn load_json(path: impl Into<PathBuf>) -> Result<Self, InteractionMatrixError> {
        let path = path.into();
        let bytes = fs::read(&path).map_err(|source| InteractionMatrixError::Io {
            path: path.clone(),
            source,
        })?;
        Self::from_json_bytes(
            bytes,
            path.clone(),
            InteractionProvenance::JsonFile { path },
        )
    }

    pub fn from_generated(
        values: DenseMatrix<f64>,
        generator: GeneratorProvenance,
    ) -> Result<Self, InteractionMatrixError> {
        Self::resolve(
            Arc::new(values),
            InteractionProvenance::Generated { generator },
        )
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
    /// Apply this interaction matrix to a contiguous batch of species vectors.
    #[inline]
    pub fn mul_vectors_into(&self, input: &[f64], output: &mut [f64]) -> Result<(), MatrixError> {
        self.values.mul_vectors_into(input, output)
    }
    pub const fn provenance(&self) -> &InteractionProvenance {
        &self.provenance
    }
    /// Return `A - A^T` from this exact matrix realization.
    pub fn antisymmetrize(&self) -> Result<Self, InteractionMatrixError> {
        let transposed = self.values.transpose();
        let values = self.values.sub(&transposed);
        self.derive(values, InteractionTransformation::Antisymmetrize)
    }

    /// Return this matrix multiplied by a finite scalar.
    pub fn scale(&self, scalar: f64) -> Result<Self, InteractionMatrixError> {
        require_finite_transform("scalar", scalar)?;
        self.derive(
            self.values.scalar_mul(scalar),
            InteractionTransformation::Scale { scalar },
        )
    }

    /// Return the elementwise absolute value of this matrix.
    pub fn abs(&self) -> Result<Self, InteractionMatrixError> {
        self.derive(self.values.abs(), InteractionTransformation::Abs)
    }

    /// Raise every entry below `minimum` to that finite lower bound.
    pub fn clamp_min(&self, minimum: f64) -> Result<Self, InteractionMatrixError> {
        require_finite_transform("minimum", minimum)?;
        let rows = self.values.rows();
        let columns = self.values.cols();
        self.derive(
            DenseMatrix::from_fn(rows, columns, |row, column| {
                self.values.get(row as isize, column as isize).max(minimum)
            }),
            InteractionTransformation::ClampMin { minimum },
        )
    }

    /// Lower every entry above `maximum` to that finite upper bound.
    pub fn clamp_max(&self, maximum: f64) -> Result<Self, InteractionMatrixError> {
        require_finite_transform("maximum", maximum)?;
        let rows = self.values.rows();
        let columns = self.values.cols();
        self.derive(
            DenseMatrix::from_fn(rows, columns, |row, column| {
                self.values.get(row as isize, column as isize).min(maximum)
            }),
            InteractionTransformation::ClampMax { maximum },
        )
    }

    /// Require every absolute entry to be at most `threshold`.
    ///
    /// The maximum is reduced by PiP's backend-native parallel implementation.
    pub fn ensure_max_abs_at_most(&self, threshold: f64) -> Result<(), InteractionMatrixError> {
        if !threshold.is_finite() || threshold < 0.0 {
            return Err(InteractionMatrixError::InvalidTransformationParameter {
                name: "threshold",
                value: threshold,
            });
        }
        let maximum = self.values.max_abs_real();
        if maximum > threshold {
            return Err(InteractionMatrixError::MaximumAbsoluteEntryExceeded {
                threshold,
                maximum,
            });
        }
        Ok(())
    }

    /// Scale down uniformly until every absolute entry is at most `threshold`.
    ///
    /// A matrix already within the threshold is not enlarged.
    pub fn normalize(&self, threshold: f64) -> Result<Self, InteractionMatrixError> {
        if !threshold.is_finite() || threshold < 0.0 {
            return Err(InteractionMatrixError::InvalidTransformationParameter {
                name: "threshold",
                value: threshold,
            });
        }
        let maximum = self.values.max_abs_real();
        let scalar = if maximum > threshold {
            threshold / maximum
        } else {
            1.0
        };
        let transformation = InteractionTransformation::Normalize {
            threshold,
            maximum,
            scalar,
        };
        if scalar == 1.0 {
            return Ok(Self {
                values: Arc::clone(&self.values),
                provenance: self.derived_provenance(transformation),
            });
        }
        self.derive(self.values.scalar_mul(scalar), transformation)
    }
    pub fn generator_rng_record(&self) -> Result<Option<RngRecord>, RngRecordError> {
        self.provenance.generator_rng_record()
    }

    fn from_json_bytes(
        bytes: Vec<u8>,
        path: PathBuf,
        provenance: InteractionProvenance,
    ) -> Result<Self, InteractionMatrixError> {
        let values =
            serde_json::from_slice(&bytes).map_err(|source| InteractionMatrixError::Json {
                path: path.clone(),
                source,
            })?;
        Self::resolve(Arc::new(values), provenance)
    }

    fn derive(
        &self,
        values: DenseMatrix<f64>,
        transformation: InteractionTransformation,
    ) -> Result<Self, InteractionMatrixError> {
        Self::resolve(Arc::new(values), self.derived_provenance(transformation))
    }

    fn derived_provenance(
        &self,
        transformation: InteractionTransformation,
    ) -> InteractionProvenance {
        InteractionProvenance::Derived {
            source: Box::new(self.provenance.clone()),
            transformation,
        }
    }

    fn resolve(
        values: Arc<DenseMatrix<f64>>,
        provenance: InteractionProvenance,
    ) -> Result<Self, InteractionMatrixError> {
        let rows = values.rows();
        let columns = values.cols();
        if rows != columns {
            return Err(InteractionMatrixError::NonSquare { rows, columns });
        }
        for flat in 0..values.size() {
            let value = values.get_flat(flat as isize);
            if !value.is_finite() {
                return Err(InteractionMatrixError::NonFiniteEntry {
                    row: flat / columns,
                    column: flat % columns,
                    value,
                });
            }
        }
        Ok(Self { values, provenance })
    }
}

fn require_finite_transform(name: &'static str, value: f64) -> Result<(), InteractionMatrixError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(InteractionMatrixError::InvalidTransformationParameter { name, value })
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum InteractionSourceKind {
    InMemory,
    Inline,
    JsonFile,
    Generated,
    Derived,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "operation", rename_all = "snake_case", deny_unknown_fields)]
pub enum InteractionTransformation {
    Antisymmetrize,
    Abs,
    Scale {
        scalar: f64,
    },
    ClampMin {
        minimum: f64,
    },
    ClampMax {
        maximum: f64,
    },
    Normalize {
        threshold: f64,
        maximum: f64,
        scalar: f64,
    },
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum InteractionProvenance {
    InMemory {
        label: Option<String>,
    },
    Inline,
    JsonFile {
        path: PathBuf,
    },
    Generated {
        generator: GeneratorProvenance,
    },
    Derived {
        source: Box<InteractionProvenance>,
        transformation: InteractionTransformation,
    },
}

impl InteractionProvenance {
    pub const fn kind(&self) -> InteractionSourceKind {
        match self {
            Self::InMemory { .. } => InteractionSourceKind::InMemory,
            Self::Inline => InteractionSourceKind::Inline,
            Self::JsonFile { .. } => InteractionSourceKind::JsonFile,
            Self::Generated { .. } => InteractionSourceKind::Generated,
            Self::Derived { .. } => InteractionSourceKind::Derived,
        }
    }
    pub const fn generator(&self) -> Option<&GeneratorProvenance> {
        match self {
            Self::Generated { generator } => Some(generator),
            Self::Derived { source, .. } => source.generator(),
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
    provenance: InteractionProvenance,
}

/// Portable pointer to a verified interaction artifact produced by an earlier execution.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct InteractionArtifactReference {
    format: String,
    execution_directory: PathBuf,
    descriptor: InteractionArtifactDescriptor,
}

impl InteractionArtifactReference {
    pub fn new(
        execution_directory: impl Into<PathBuf>,
        descriptor: InteractionArtifactDescriptor,
    ) -> Self {
        Self {
            format: INTERACTION_ARTIFACT_REFERENCE_FORMAT.to_owned(),
            execution_directory: execution_directory.into(),
            descriptor,
        }
    }

    pub fn load_json(path: impl Into<PathBuf>) -> Result<Self, InteractionArtifactLoadError> {
        let path = path.into();
        let bytes = fs::read(&path).map_err(|source| InteractionMatrixError::Io {
            path: path.clone(),
            source,
        })?;
        let reference: Self =
            serde_json::from_slice(&bytes).map_err(|source| InteractionMatrixError::Json {
                path: path.clone(),
                source,
            })?;
        if reference.format != INTERACTION_ARTIFACT_REFERENCE_FORMAT
            || reference.execution_directory.as_os_str().is_empty()
        {
            return Err(InteractionArtifactLoadError::InvalidDescriptor);
        }
        Ok(reference)
    }

    pub fn to_json_bytes(&self) -> Result<Vec<u8>, serde_json::Error> {
        serde_json::to_vec_pretty(self)
    }

    pub fn resolve(&self) -> Result<InteractionMatrix, InteractionArtifactLoadError> {
        load_verified_interaction_matrix(&self.execution_directory, &self.descriptor)
    }

    pub fn execution_directory(&self) -> &Path {
        &self.execution_directory
    }

    pub const fn descriptor(&self) -> &InteractionArtifactDescriptor {
        &self.descriptor
    }
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
        self.provenance.kind()
    }
    pub const fn generator(&self) -> Option<&GeneratorProvenance> {
        self.provenance.generator()
    }
    pub const fn provenance(&self) -> &InteractionProvenance {
        &self.provenance
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
            provenance: matrix.provenance().clone(),
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
    let matrix = InteractionMatrix::from_json_bytes(
        verified.into_bytes(),
        path,
        descriptor.provenance.clone(),
    )?;
    if matrix.species() != descriptor.species {
        return Err(InteractionMatrixError::SpeciesMismatch {
            expected: descriptor.species,
            actual: matrix.species(),
        }
        .into());
    }
    Ok(matrix)
}

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum InteractionRecipeError {
    #[error("interaction recipe requires at least one species")]
    EmptySpecies,
    #[error("interaction recipe parameter {name} is invalid: {value}")]
    InvalidParameter { name: &'static str, value: f64 },
    #[error("uniform interaction range requires minimum < maximum, got [{minimum}, {maximum})")]
    InvalidRange { minimum: f64, maximum: f64 },
    #[error("interaction family {family} does not support sampled diagonal entries")]
    SampledDiagonalUnsupported { family: &'static str },
    #[error(transparent)]
    Rng(#[from] RngConfigError),
    #[error(transparent)]
    TensorRand(#[from] TensorRandError),
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
    #[error("interaction matrix transformation parameter {name} is invalid: {value}")]
    InvalidTransformationParameter { name: &'static str, value: f64 },
    #[error("interaction matrix maximum absolute entry {maximum} exceeds threshold {threshold}")]
    MaximumAbsoluteEntryExceeded { threshold: f64, maximum: f64 },
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
