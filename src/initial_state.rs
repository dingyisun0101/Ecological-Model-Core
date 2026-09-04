//! Reproducible categorical ecological initial states and verified artifacts.

use std::fs;
use std::path::{Path, PathBuf};

use physics_in_parallel::prelude::basic::{
    ResolvedRng, RngMethod, SquareLattice, SquareLatticeGeometry, SquareLatticeGeometryError,
    SquareLatticeInitMethod,
};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use thiserror::Error;

use crate::artifact::{
    ArtifactDescriptor, ArtifactDisposition, ArtifactError, ArtifactLoadError,
    load_verified_artifact, persist_artifact,
};

pub const INITIAL_STATE_FORMAT: &str = "ecological.initial-state.v2";
pub const INITIAL_STATE_METADATA_KEY: &str = "initial_state";
pub const INITIAL_STATE_ARTIFACT_REFERENCE_FORMAT: &str =
    "ecological.initial-state-artifact-reference.v2";

pub type CategoricalSpace = SquareLattice<usize>;
pub type TaxonCounts = Vec<usize>;

/// Source of relative taxon weights for categorical placement.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum DistributionSource {
    Uniform,
    Inline { weights: Vec<f64> },
    Json { path: PathBuf },
}

impl DistributionSource {
    pub fn path(&self) -> Option<&Path> {
        match self {
            Self::Json { path } => Some(path),
            Self::Uniform | Self::Inline { .. } => None,
        }
    }

    fn validate(&self, num_taxa: usize) -> Result<(), InitialStateError> {
        match self {
            Self::Uniform => Ok(()),
            Self::Inline { weights } => validate_distribution(weights, num_taxa),
            Self::Json { path } if path.as_os_str().is_empty() => Err(InitialStateError::EmptyPath),
            Self::Json { .. } => Ok(()),
        }
    }

    fn resolve(self, num_taxa: usize) -> Result<Option<Vec<f64>>, InitialStateError> {
        let weights = match self {
            Self::Uniform => return Ok(None),
            Self::Inline { weights } => weights,
            Self::Json { path } => {
                let bytes = fs::read(&path).map_err(|source| InitialStateError::Io {
                    operation: "read distribution",
                    path: path.clone(),
                    source,
                })?;
                serde_json::from_slice(&bytes)
                    .map_err(|source| InitialStateError::Json { path, source })?
            }
        };
        validate_distribution(&weights, num_taxa)?;
        Ok(Some(weights))
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum InitializationMethod {
    Random,
    BalancedUniform,
    CenteredSeed,
    CenteredDominantSeed,
}

/// Reproducible scientific instructions for one categorical state.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "method", rename_all = "snake_case", deny_unknown_fields)]
pub enum InitialStateRecipe {
    Random {
        distribution: DistributionSource,
        #[serde(default = "default_indexed_rng")]
        rng: ResolvedRng,
    },
    /// Assign the minimum-variance uniform taxon counts, then shuffle all sites.
    BalancedUniform {
        #[serde(default = "default_indexed_rng")]
        rng: ResolvedRng,
    },
    CenteredSeed {
        distribution: DistributionSource,
        seed_taxon: usize,
        seed_radius: usize,
        #[serde(default = "default_indexed_rng")]
        rng: ResolvedRng,
    },
    CenteredDominantSeed {
        distribution: DistributionSource,
        seed_radius: usize,
        #[serde(default = "default_indexed_rng")]
        rng: ResolvedRng,
    },
}

impl InitialStateRecipe {
    pub const fn method(&self) -> InitializationMethod {
        match self {
            Self::Random { .. } => InitializationMethod::Random,
            Self::BalancedUniform { .. } => InitializationMethod::BalancedUniform,
            Self::CenteredSeed { .. } => InitializationMethod::CenteredSeed,
            Self::CenteredDominantSeed { .. } => InitializationMethod::CenteredDominantSeed,
        }
    }

    pub fn distribution_path(&self) -> Option<&Path> {
        match self {
            Self::Random { distribution, .. }
            | Self::CenteredSeed { distribution, .. }
            | Self::CenteredDominantSeed { distribution, .. } => distribution.path(),
            Self::BalancedUniform { .. } => None,
        }
    }

    pub fn validate(
        &self,
        lattice: &SquareLatticeGeometry,
        num_taxa: usize,
    ) -> Result<(), InitialStateError> {
        if num_taxa == 0 {
            return Err(InitialStateError::ZeroTaxa);
        }
        match self {
            Self::Random { distribution, .. } => distribution.validate(num_taxa),
            Self::BalancedUniform { .. } => Ok(()),
            Self::CenteredSeed {
                distribution,
                seed_taxon,
                seed_radius,
                ..
            } => {
                distribution.validate(num_taxa)?;
                validate_seed(lattice, num_taxa, *seed_taxon, *seed_radius)
            }
            Self::CenteredDominantSeed {
                distribution,
                seed_radius,
                ..
            } => {
                distribution.validate(num_taxa)?;
                validate_seed_geometry(lattice, *seed_radius)
            }
        }
    }

    pub fn create(
        self,
        lattice: SquareLatticeGeometry,
        num_taxa: usize,
    ) -> Result<InitialState, InitialStateError> {
        self.validate(&lattice, num_taxa)?;
        let method = self.method();
        let (distribution, rng, explicit_seed, seed_radius, dominant) = match self {
            Self::Random { distribution, rng } => (distribution, rng, None, None, false),
            Self::CenteredSeed {
                distribution,
                seed_taxon,
                seed_radius,
                rng,
            } => (
                distribution,
                rng,
                Some(seed_taxon),
                Some(seed_radius),
                false,
            ),
            Self::CenteredDominantSeed {
                distribution,
                seed_radius,
                rng,
            } => (distribution, rng, None, Some(seed_radius), true),
            Self::BalancedUniform { rng } => {
                let values = balanced_uniform_values(lattice.num_sites(), num_taxa);
                let space = CategoricalSpace::new(
                    lattice,
                    SquareLatticeInitMethod::ShuffledValues { values, rng },
                )?;
                let counts = count_taxa(&space, num_taxa)?;
                let rng = Some(resolved_rng_from_space(&space)?);
                return Ok(InitialState {
                    num_taxa,
                    method,
                    seed_taxon: None,
                    rng,
                    space,
                    counts,
                });
            }
        };
        let mut weights = distribution.resolve(num_taxa)?;
        let seed_taxon = if dominant {
            let background = weights.get_or_insert_with(|| vec![1.0; num_taxa]);
            Some(remove_dominant_taxon(background)?)
        } else {
            explicit_seed
        };
        let mut space = CategoricalSpace::new(
            lattice,
            SquareLatticeInitMethod::RandomChoices {
                choices: (0..num_taxa).collect(),
                weights,
                rng,
            },
        )?;
        let mut counts = count_taxa(&space, num_taxa)?;
        if let (Some(taxon), Some(radius)) = (seed_taxon, seed_radius) {
            let shape = space.geometry().shape().to_vec();
            plant_centered_seed(&mut space, &mut counts, &shape, taxon, radius);
        }
        let rng = Some(resolved_rng_from_space(&space)?);
        Ok(InitialState {
            num_taxa,
            method,
            seed_taxon,
            rng,
            space,
            counts,
        })
    }
}

/// A generated recipe or a verified prior-execution artifact.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "source", rename_all = "snake_case", deny_unknown_fields)]
pub enum InitialStateSource {
    Recipe {
        recipe: InitialStateRecipe,
    },
    VerifiedArtifact {
        artifact_root: PathBuf,
        descriptor: InitialStateArtifactDescriptor,
    },
}

impl InitialStateSource {
    pub fn validate(
        &self,
        lattice: &SquareLatticeGeometry,
        num_taxa: usize,
    ) -> Result<(), InitialStateError> {
        match self {
            Self::Recipe { recipe } => recipe.validate(lattice, num_taxa),
            Self::VerifiedArtifact {
                artifact_root,
                descriptor,
            } => {
                if artifact_root.as_os_str().is_empty() {
                    return Err(InitialStateError::EmptyPath);
                }
                if descriptor.lattice() != lattice {
                    return Err(InitialStateError::LatticeMismatch);
                }
                if descriptor.num_taxa() != num_taxa {
                    return Err(InitialStateError::TaxonDimensionMismatch {
                        expected: num_taxa,
                        actual: descriptor.num_taxa(),
                    });
                }
                Ok(())
            }
        }
    }

    pub fn resolve(
        &self,
        lattice: SquareLatticeGeometry,
        num_taxa: usize,
    ) -> Result<InitialState, InitialStateError> {
        self.validate(&lattice, num_taxa)?;
        match self {
            Self::Recipe { recipe } => recipe.clone().create(lattice, num_taxa),
            Self::VerifiedArtifact {
                artifact_root,
                descriptor,
            } => load_verified_initial_state(artifact_root, descriptor),
        }
    }
}

#[derive(Debug)]
pub struct InitialState {
    num_taxa: usize,
    method: InitializationMethod,
    seed_taxon: Option<usize>,
    rng: Option<ResolvedRng>,
    space: CategoricalSpace,
    counts: TaxonCounts,
}

impl InitialState {
    pub const fn num_taxa(&self) -> usize {
        self.num_taxa
    }
    pub const fn method(&self) -> InitializationMethod {
        self.method
    }
    pub const fn seed_taxon(&self) -> Option<usize> {
        self.seed_taxon
    }
    /// Returns the fully resolved PiP randomness used to construct the state.
    pub const fn rng_config(&self) -> Option<ResolvedRng> {
        self.rng
    }
    pub const fn space(&self) -> &CategoricalSpace {
        &self.space
    }
    pub fn counts(&self) -> &[usize] {
        &self.counts
    }
    /// Returns the exact aggregate relative frequencies represented by the lattice.
    pub fn frequencies(&self) -> Vec<f64> {
        let total = self.space.num_sites() as f64;
        self.counts
            .iter()
            .map(|&count| count as f64 / total)
            .collect()
    }
    /// Encodes the complete reproducible state in eco_core's canonical JSON format.
    pub fn to_json_bytes(&self) -> Result<Vec<u8>, serde_json::Error> {
        serde_json::to_vec(&InitialStateDocumentRef::from(self))
    }
    pub fn clone_space(&self) -> CategoricalSpace {
        self.space.clone()
    }
    pub fn into_parts(self) -> (CategoricalSpace, TaxonCounts) {
        (self.space, self.counts)
    }
}

#[derive(Serialize)]
struct InitialStateDocumentRef<'a> {
    format: &'static str,
    num_taxa: usize,
    method: InitializationMethod,
    seed_taxon: Option<usize>,
    rng: Option<ResolvedRng>,
    lattice: &'a SquareLatticeGeometry,
    sites: &'a [usize],
}

impl<'a> From<&'a InitialState> for InitialStateDocumentRef<'a> {
    fn from(initial: &'a InitialState) -> Self {
        Self {
            format: INITIAL_STATE_FORMAT,
            num_taxa: initial.num_taxa,
            method: initial.method,
            seed_taxon: initial.seed_taxon,
            rng: initial.rng,
            lattice: initial.space.geometry(),
            sites: initial.space.data(),
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct InitialStateArtifactDescriptor {
    format: String,
    num_taxa: usize,
    lattice: SquareLatticeGeometry,
    method: InitializationMethod,
    #[serde(skip_serializing_if = "Option::is_none")]
    seed_taxon: Option<usize>,
    #[serde(flatten)]
    artifact: ArtifactDescriptor,
}

/// Serializable pointer to a verified initial-state artifact beneath an authored root.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct InitialStateArtifactReference {
    format: String,
    artifact_root: PathBuf,
    descriptor: InitialStateArtifactDescriptor,
}

impl InitialStateArtifactReference {
    pub fn new(artifact_root: PathBuf, descriptor: InitialStateArtifactDescriptor) -> Self {
        Self {
            format: INITIAL_STATE_ARTIFACT_REFERENCE_FORMAT.to_owned(),
            artifact_root,
            descriptor,
        }
    }

    pub fn from_json_bytes(bytes: &[u8]) -> Result<Self, InitialStateError> {
        let reference: Self =
            serde_json::from_slice(bytes).map_err(|source| InitialStateError::Json {
                path: PathBuf::from("<initial-state-artifact-reference>"),
                source,
            })?;
        reference.validate()?;
        Ok(reference)
    }

    pub fn load_json(path: &Path) -> Result<Self, InitialStateError> {
        let path = path.to_path_buf();
        let bytes = fs::read(&path).map_err(|source| InitialStateError::Io {
            operation: "read initial-state artifact reference",
            path: path.clone(),
            source,
        })?;
        let reference: Self =
            serde_json::from_slice(&bytes).map_err(|source| InitialStateError::Json {
                path: path.clone(),
                source,
            })?;
        reference.validate()?;
        Ok(reference)
    }

    pub fn to_json_bytes(&self) -> Result<Vec<u8>, serde_json::Error> {
        serde_json::to_vec_pretty(self)
    }

    pub fn resolve(&self) -> Result<InitialState, InitialStateError> {
        self.validate()?;
        load_verified_initial_state(&self.artifact_root, &self.descriptor)
    }

    /// Validates the reference envelope without reading the artifact.
    pub fn validate(&self) -> Result<(), InitialStateError> {
        if self.format != INITIAL_STATE_ARTIFACT_REFERENCE_FORMAT
            || self.artifact_root.as_os_str().is_empty()
        {
            Err(InitialStateError::InvalidArtifactReference)
        } else {
            Ok(())
        }
    }

    pub fn artifact_root(&self) -> &Path {
        &self.artifact_root
    }

    pub const fn descriptor(&self) -> &InitialStateArtifactDescriptor {
        &self.descriptor
    }
}

impl InitialStateArtifactDescriptor {
    pub fn format(&self) -> &str {
        &self.format
    }
    pub const fn num_taxa(&self) -> usize {
        self.num_taxa
    }
    pub const fn lattice(&self) -> &SquareLatticeGeometry {
        &self.lattice
    }
    pub const fn method(&self) -> InitializationMethod {
        self.method
    }
    pub const fn seed_taxon(&self) -> Option<usize> {
        self.seed_taxon
    }
    pub fn sha256(&self) -> &str {
        self.artifact.sha256()
    }
    pub fn path(&self) -> &Path {
        self.artifact.path()
    }
    pub fn insert_into_metadata(&self, metadata: &mut Map<String, Value>) -> Option<Value> {
        metadata.insert(
            INITIAL_STATE_METADATA_KEY.to_owned(),
            serde_json::to_value(self).expect("initial-state descriptor is JSON-compatible"),
        )
    }
}

#[derive(Clone, Debug)]
pub struct PersistedInitialState {
    descriptor: InitialStateArtifactDescriptor,
    disposition: ArtifactDisposition,
}

impl PersistedInitialState {
    pub const fn descriptor(&self) -> &InitialStateArtifactDescriptor {
        &self.descriptor
    }
    pub const fn disposition(&self) -> ArtifactDisposition {
        self.disposition
    }
    pub fn into_descriptor(self) -> InitialStateArtifactDescriptor {
        self.descriptor
    }
}

pub fn persist_initial_state(
    artifact_root: &Path,
    initial: &InitialState,
) -> Result<PersistedInitialState, InitialStateError> {
    let bytes = serde_json::to_vec(&InitialStateDocumentRef::from(initial))?;
    let persisted = persist_artifact(artifact_root, "initial-state", "json", &bytes)?;
    Ok(PersistedInitialState {
        descriptor: InitialStateArtifactDescriptor {
            format: INITIAL_STATE_FORMAT.to_owned(),
            num_taxa: initial.num_taxa,
            lattice: initial.space.geometry().clone(),
            method: initial.method,
            seed_taxon: initial.seed_taxon,
            artifact: persisted.descriptor().clone(),
        },
        disposition: persisted.disposition(),
    })
}

pub fn load_verified_initial_state(
    artifact_root: &Path,
    descriptor: &InitialStateArtifactDescriptor,
) -> Result<InitialState, InitialStateError> {
    if descriptor.format != INITIAL_STATE_FORMAT {
        return Err(InitialStateError::UnsupportedFormat {
            actual: descriptor.format.clone(),
        });
    }
    let verified = load_verified_artifact(artifact_root, &descriptor.artifact)?;
    let document: InitialStateDocument =
        serde_json::from_slice(verified.bytes()).map_err(|source| InitialStateError::Json {
            path: verified.path().to_path_buf(),
            source,
        })?;
    let initial = document.resolve(descriptor.lattice.clone(), descriptor.num_taxa)?;
    if initial.method != descriptor.method || initial.seed_taxon != descriptor.seed_taxon {
        return Err(InitialStateError::DescriptorMismatch);
    }
    Ok(initial)
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct InitialStateDocument {
    format: String,
    num_taxa: usize,
    method: InitializationMethod,
    seed_taxon: Option<usize>,
    rng: Option<ResolvedRng>,
    lattice: SquareLatticeGeometry,
    sites: Vec<usize>,
}

impl InitialStateDocument {
    fn resolve(
        self,
        expected_lattice: SquareLatticeGeometry,
        expected_num_taxa: usize,
    ) -> Result<InitialState, InitialStateError> {
        if self.format != INITIAL_STATE_FORMAT {
            return Err(InitialStateError::UnsupportedFormat {
                actual: self.format,
            });
        }
        if self.num_taxa != expected_num_taxa {
            return Err(InitialStateError::TaxonDimensionMismatch {
                expected: expected_num_taxa,
                actual: self.num_taxa,
            });
        }
        if self.lattice != expected_lattice {
            return Err(InitialStateError::LatticeMismatch);
        }
        if self.seed_taxon.is_some_and(|taxon| taxon >= self.num_taxa) {
            return Err(InitialStateError::SeedTaxonOutOfRange {
                seed_taxon: self.seed_taxon.expect("checked Some"),
                num_taxa: self.num_taxa,
            });
        }
        let space = CategoricalSpace::new(
            self.lattice,
            SquareLatticeInitMethod::Values { values: self.sites },
        )?;
        let counts = count_taxa(&space, self.num_taxa)?;
        Ok(InitialState {
            num_taxa: self.num_taxa,
            method: self.method,
            seed_taxon: self.seed_taxon,
            rng: self.rng,
            space,
            counts,
        })
    }
}

fn resolved_rng_from_space(space: &CategoricalSpace) -> Result<ResolvedRng, InitialStateError> {
    space
        .initialization_resolved_rng()
        .ok_or(InitialStateError::MissingResolvedRngConfig)
}

fn default_indexed_rng() -> ResolvedRng {
    ResolvedRng::from_entropy(RngMethod::IndexedSplitMix64)
}

fn validate_distribution(weights: &[f64], num_taxa: usize) -> Result<(), InitialStateError> {
    if weights.len() != num_taxa {
        return Err(InitialStateError::DistributionLength {
            expected: num_taxa,
            actual: weights.len(),
        });
    }
    let mut total = 0.0;
    for (taxon, &weight) in weights.iter().enumerate() {
        if !weight.is_finite() || weight < 0.0 {
            return Err(InitialStateError::InvalidWeight { taxon, weight });
        }
        total += weight;
    }
    if !total.is_finite() || total <= 0.0 {
        return Err(InitialStateError::NonPositiveDistribution);
    }
    Ok(())
}

fn remove_dominant_taxon(weights: &mut [f64]) -> Result<usize, InitialStateError> {
    validate_distribution(weights, weights.len())?;
    let dominant = weights
        .iter()
        .enumerate()
        .max_by(|left, right| left.1.total_cmp(right.1).then(right.0.cmp(&left.0)))
        .map(|(taxon, _)| taxon)
        .expect("validated distribution is nonempty");
    let background_total = weights
        .iter()
        .enumerate()
        .filter_map(|(taxon, weight)| (taxon != dominant).then_some(*weight))
        .sum::<f64>();
    if background_total <= 0.0 {
        return Err(InitialStateError::NoDominantSeedBackground);
    }
    weights[dominant] = 0.0;
    for weight in weights {
        *weight /= background_total;
    }
    Ok(dominant)
}

fn count_taxa(space: &CategoricalSpace, num_taxa: usize) -> Result<TaxonCounts, InitialStateError> {
    let mut counts = vec![0; num_taxa];
    for (site, &taxon) in space.data().iter().enumerate() {
        let count = counts
            .get_mut(taxon)
            .ok_or(InitialStateError::SpaceTaxonOutOfRange {
                site,
                taxon,
                num_taxa,
            })?;
        *count += 1;
    }
    Ok(counts)
}

fn balanced_uniform_values(num_sites: usize, num_taxa: usize) -> Vec<usize> {
    let per_taxon = num_sites / num_taxa;
    let remainder = num_sites % num_taxa;
    let mut values = Vec::with_capacity(num_sites);
    for taxon in 0..num_taxa {
        values.extend(std::iter::repeat_n(
            taxon,
            per_taxon + usize::from(taxon < remainder),
        ));
    }
    debug_assert_eq!(values.len(), num_sites);
    values
}

fn validate_seed(
    lattice: &SquareLatticeGeometry,
    num_taxa: usize,
    seed_taxon: usize,
    seed_radius: usize,
) -> Result<(), InitialStateError> {
    if seed_taxon >= num_taxa {
        return Err(InitialStateError::SeedTaxonOutOfRange {
            seed_taxon,
            num_taxa,
        });
    }
    validate_seed_geometry(lattice, seed_radius)
}

fn validate_seed_geometry(
    lattice: &SquareLatticeGeometry,
    seed_radius: usize,
) -> Result<(), InitialStateError> {
    let width = seed_radius
        .checked_mul(2)
        .and_then(|diameter| diameter.checked_add(1))
        .ok_or(InitialStateError::SeedWidthOverflow {
            radius: seed_radius,
        })?;
    for (axis, &length) in lattice.shape().iter().enumerate() {
        if width > length {
            return Err(InitialStateError::SeedDoesNotFit {
                axis,
                width,
                length,
            });
        }
    }
    Ok(())
}

fn plant_centered_seed(
    space: &mut CategoricalSpace,
    counts: &mut [usize],
    shape: &[usize],
    seed_taxon: usize,
    radius: usize,
) {
    let starts = shape
        .iter()
        .map(|length| length / 2 - radius)
        .collect::<Vec<_>>();
    let widths = vec![radius * 2 + 1; shape.len()];
    let seed_sites = widths.iter().product::<usize>();
    let mut coordinate = vec![0isize; shape.len()];
    for local_flat in 0..seed_sites {
        let mut local = local_flat;
        for axis in (0..shape.len()).rev() {
            coordinate[axis] = (starts[axis] + local % widths[axis]) as isize;
            local /= widths[axis];
        }
        let previous = *space.get(&coordinate);
        if previous != seed_taxon {
            counts[previous] -= 1;
            counts[seed_taxon] += 1;
            space.set(&coordinate, seed_taxon);
        }
    }
}

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum InitialStateError {
    #[error("num_taxa must be positive")]
    ZeroTaxa,
    #[error("distribution length is {actual}, expected {expected}")]
    DistributionLength { expected: usize, actual: usize },
    #[error("distribution weight {taxon} is invalid: {weight}")]
    InvalidWeight { taxon: usize, weight: f64 },
    #[error("distribution must have positive finite total weight")]
    NonPositiveDistribution,
    #[error("dominant-seed initialization requires positive background mass")]
    NoDominantSeedBackground,
    #[error("seed taxon {seed_taxon} is outside 0..{num_taxa}")]
    SeedTaxonOutOfRange { seed_taxon: usize, num_taxa: usize },
    #[error("seed radius {radius} is too large")]
    SeedWidthOverflow { radius: usize },
    #[error("seed width {width} does not fit lattice axis {axis} of length {length}")]
    SeedDoesNotFit {
        axis: usize,
        width: usize,
        length: usize,
    },
    #[error("path must not be empty")]
    EmptyPath,
    #[error("initial-state format `{actual}` is unsupported")]
    UnsupportedFormat { actual: String },
    #[error("initial state declares {actual} taxa, expected {expected}")]
    TaxonDimensionMismatch { expected: usize, actual: usize },
    #[error("initial-state lattice does not match the expected lattice")]
    LatticeMismatch,
    #[error("initial-state descriptor does not match its verified document")]
    DescriptorMismatch,
    #[error("invalid initial-state artifact reference")]
    InvalidArtifactReference,
    #[error("initial-state site {site} contains taxon {taxon}, outside 0..{num_taxa}")]
    SpaceTaxonOutOfRange {
        site: usize,
        taxon: usize,
        num_taxa: usize,
    },
    #[error("failed to {operation} at `{path}`")]
    Io {
        operation: &'static str,
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("invalid JSON initial state at `{path}`")]
    Json {
        path: PathBuf,
        #[source]
        source: serde_json::Error,
    },
    #[error(transparent)]
    Serialize(#[from] serde_json::Error),
    #[error(transparent)]
    Lattice(#[from] SquareLatticeGeometryError),
    #[error("initialized lattice has no resolved RNG configuration")]
    MissingResolvedRngConfig,
    #[error(transparent)]
    Artifact(#[from] ArtifactError),
    #[error(transparent)]
    ArtifactLoad(#[from] ArtifactLoadError),
}
