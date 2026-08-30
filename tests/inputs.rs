use ecological_state_toolkit::initial_state::{
    InitialStateArtifactReference, InitialStateRecipe, persist_initial_state,
};
use ecological_state_toolkit::inputs::{EcologicalInputs, EcologicalInputsError};
use ecological_state_toolkit::interaction::{
    InteractionArtifactReference, InteractionMatrix, persist_interaction_matrix,
};
use physics_in_parallel::prelude::basic::{RngConfig, SquareLatticeConfig};

fn prepared_inputs(
    root: &std::path::Path,
    interaction_species: usize,
    initial_taxa: usize,
) -> (InteractionArtifactReference, InitialStateArtifactReference) {
    let matrix =
        InteractionMatrix::from_rows(vec![vec![0.0; interaction_species]; interaction_species])
            .unwrap();
    let interaction = persist_interaction_matrix(root, &matrix).unwrap();
    let initial = InitialStateRecipe::BalancedUniform {
        rng: RngConfig::new(Some(81), None),
    }
    .create(SquareLatticeConfig::periodic(&[12]), initial_taxa)
    .unwrap();
    let initial = persist_initial_state(root, &initial).unwrap();
    (
        InteractionArtifactReference::new(root.to_path_buf(), interaction.into_descriptor()),
        InitialStateArtifactReference::new(root.to_path_buf(), initial.into_descriptor()),
    )
}

#[test]
fn shared_envelope_round_trips_and_resolves_both_verified_inputs() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path().join("prepared");
    let (interaction, initial_state) = prepared_inputs(&root, 3, 3);
    let inputs = EcologicalInputs::new(interaction, initial_state).unwrap();

    let encoded = serde_json::to_vec(&inputs).unwrap();
    let decoded: EcologicalInputs = serde_json::from_slice(&encoded).unwrap();
    assert_eq!(decoded, inputs);
    let resolved = decoded.resolve().unwrap();
    assert_eq!(resolved.interaction().species(), 3);
    assert_eq!(resolved.initial_state().num_taxa(), 3);
    assert_eq!(resolved.initial_state().counts().iter().sum::<usize>(), 12);
    assert_eq!(
        resolved.initial_state().frequencies().iter().sum::<f64>(),
        1.0
    );
}

#[test]
fn shared_envelope_rejects_dimension_mismatch_before_artifact_io() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path().join("prepared");
    let (interaction, initial_state) = prepared_inputs(&root, 2, 3);
    assert!(matches!(
        EcologicalInputs::new(interaction, initial_state),
        Err(EcologicalInputsError::TaxonDimensionMismatch {
            interaction_species: 2,
            initial_taxa: 3,
        })
    ));
}

#[test]
fn deserialization_rejects_an_invalid_reference_envelope() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path().join("prepared");
    let (interaction, initial_state) = prepared_inputs(&root, 2, 2);
    let inputs = EcologicalInputs::new(interaction, initial_state).unwrap();
    let mut document = serde_json::to_value(inputs).unwrap();
    document["initial_state"]["artifact_root"] = "".into();
    assert!(serde_json::from_value::<EcologicalInputs>(document).is_err());
}
