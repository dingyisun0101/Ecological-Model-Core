//! Model-neutral references to one prepared ecological realization.

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::initial_state::{InitialState, InitialStateArtifactReference, InitialStateError};
use crate::interaction::{
    InteractionArtifactLoadError, InteractionArtifactReference, InteractionMatrix,
};

/// Immutable prepared inputs consumed by one ecological model member.
///
/// Application-facing model configuration uses this envelope instead of model
/// recipes, inline matrices, independent frequency vectors, or unverified
/// paths. The interaction is already in the receiving model's convention. The
/// initial-state reference is shareable across models: a lattice model consumes
/// its categorical space and counts, while a continuous model derives exact
/// frequencies from the same artifact.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct EcologicalInputs {
    interaction: InteractionArtifactReference,
    initial_state: InitialStateArtifactReference,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct EcologicalInputsDocument {
    interaction: InteractionArtifactReference,
    initial_state: InitialStateArtifactReference,
}

impl<'de> Deserialize<'de> for EcologicalInputs {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let document = EcologicalInputsDocument::deserialize(deserializer)?;
        Self::new(document.interaction, document.initial_state).map_err(serde::de::Error::custom)
    }
}

impl EcologicalInputs {
    /// Builds one validated input envelope without reading artifact bytes.
    pub fn new(
        interaction: InteractionArtifactReference,
        initial_state: InitialStateArtifactReference,
    ) -> Result<Self, EcologicalInputsError> {
        let inputs = Self {
            interaction,
            initial_state,
        };
        inputs.validate()?;
        Ok(inputs)
    }

    /// Validates reference envelopes and their shared taxon dimension.
    ///
    /// This operation performs no filesystem IO. Full digest and semantic
    /// verification occurs in [`Self::resolve`].
    pub fn validate(&self) -> Result<(), EcologicalInputsError> {
        self.interaction.validate()?;
        self.initial_state.validate()?;
        let interaction_species = self.interaction.descriptor().species();
        let initial_taxa = self.initial_state.descriptor().num_taxa();
        if interaction_species != initial_taxa {
            return Err(EcologicalInputsError::TaxonDimensionMismatch {
                interaction_species,
                initial_taxa,
            });
        }
        Ok(())
    }

    /// Borrows the model-ready interaction reference.
    pub const fn interaction(&self) -> &InteractionArtifactReference {
        &self.interaction
    }

    /// Borrows the canonical initial-state reference.
    pub const fn initial_state(&self) -> &InitialStateArtifactReference {
        &self.initial_state
    }

    /// Loads and verifies both immutable artifacts.
    pub fn resolve(&self) -> Result<ResolvedEcologicalInputs, EcologicalInputsError> {
        self.validate()?;
        let interaction = self.interaction.resolve()?;
        let initial_state = self.initial_state.resolve()?;
        if interaction.species() != initial_state.num_taxa() {
            return Err(EcologicalInputsError::TaxonDimensionMismatch {
                interaction_species: interaction.species(),
                initial_taxa: initial_state.num_taxa(),
            });
        }
        Ok(ResolvedEcologicalInputs {
            interaction,
            initial_state,
        })
    }
}

/// Fully verified in-memory inputs for one model member.
#[derive(Debug)]
pub struct ResolvedEcologicalInputs {
    interaction: InteractionMatrix,
    initial_state: InitialState,
}

impl ResolvedEcologicalInputs {
    /// Borrows the verified model-ready interaction matrix.
    pub const fn interaction(&self) -> &InteractionMatrix {
        &self.interaction
    }

    /// Borrows the verified canonical initial state.
    pub const fn initial_state(&self) -> &InitialState {
        &self.initial_state
    }

    /// Separates the resolved values for ownership by a model.
    pub fn into_parts(self) -> (InteractionMatrix, InitialState) {
        (self.interaction, self.initial_state)
    }
}

/// Failure to validate or resolve a prepared ecological input envelope.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum EcologicalInputsError {
    /// The interaction reference or artifact is invalid.
    #[error(transparent)]
    Interaction(#[from] InteractionArtifactLoadError),
    /// The initial-state reference or artifact is invalid.
    #[error(transparent)]
    InitialState(#[from] InitialStateError),
    /// Interaction and initial state describe different ecological dimensions.
    #[error(
        "interaction has {interaction_species} species but initial state has {initial_taxa} taxa"
    )]
    TaxonDimensionMismatch {
        /// Species inferred from the interaction descriptor or matrix.
        interaction_species: usize,
        /// Taxa declared by the initial-state descriptor or state.
        initial_taxa: usize,
    },
}
