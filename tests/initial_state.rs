use std::fs;

use ecological_state_toolkit::artifact::{ArtifactDisposition, ArtifactLoadError};
use ecological_state_toolkit::initial_state::{
    DistributionSource, INITIAL_STATE_FORMAT, InitialStateArtifactReference, InitialStateError,
    InitialStateRecipe, InitializationMethod, load_verified_initial_state, persist_initial_state,
};
use physics_in_parallel::prelude::basic::{ResolvedRng, RngMethod, SquareLatticeGeometry};

fn indexed(seed: u64) -> ResolvedRng {
    ResolvedRng::new(seed, RngMethod::IndexedSplitMix64)
}

#[test]
fn recipe_is_reproducible_and_counts_are_exact() {
    let lattice = SquareLatticeGeometry::periodic(&[9, 9]).unwrap();
    let recipe = InitialStateRecipe::CenteredSeed {
        distribution: DistributionSource::Inline {
            weights: vec![0.7, 0.3],
        },
        seed_taxon: 1,
        seed_radius: 1,
        rng: indexed(42),
    };
    let first = recipe.clone().create(lattice.clone(), 2).unwrap();
    let second = recipe.create(lattice, 2).unwrap();
    assert_eq!(first.space().data(), second.space().data());
    assert_eq!(first.counts().iter().sum::<usize>(), 81);
    assert_eq!(first.method(), InitializationMethod::CenteredSeed);
    assert_eq!(first.seed_taxon(), Some(1));
}

#[test]
fn balanced_uniform_has_the_minimum_possible_count_spread() {
    let recipe = InitialStateRecipe::BalancedUniform { rng: indexed(1201) };
    let first = recipe
        .clone()
        .create(SquareLatticeGeometry::periodic(&[5, 7]).unwrap(), 8)
        .unwrap();
    let second = recipe
        .create(SquareLatticeGeometry::periodic(&[5, 7]).unwrap(), 8)
        .unwrap();
    assert_eq!(first.method(), InitializationMethod::BalancedUniform);
    assert_eq!(first.space().data(), second.space().data());
    assert_eq!(first.counts().iter().sum::<usize>(), 35);
    assert_eq!(first.counts().iter().max().unwrap(), &5);
    assert_eq!(first.counts().iter().min().unwrap(), &4);
    assert_eq!(first.frequencies().iter().sum::<f64>(), 1.0);
}

#[test]
fn dominant_recipe_uses_first_maximal_taxon() {
    let state = InitialStateRecipe::CenteredDominantSeed {
        distribution: DistributionSource::Inline {
            weights: vec![0.4, 0.4, 0.2],
        },
        seed_radius: 1,
        rng: indexed(7),
    }
    .create(SquareLatticeGeometry::periodic(&[9]).unwrap(), 3)
    .unwrap();
    assert_eq!(state.seed_taxon(), Some(0));
    assert_eq!(state.counts()[0], 3);
}

#[test]
fn dominant_recipe_retains_tiny_positive_background_mass() {
    let state = InitialStateRecipe::CenteredDominantSeed {
        distribution: DistributionSource::Inline {
            weights: vec![1.0, 1.0e-300, 2.0e-300],
        },
        seed_radius: 1,
        rng: indexed(19),
    }
    .create(SquareLatticeGeometry::periodic(&[8]).unwrap(), 3)
    .unwrap();
    assert_eq!(state.seed_taxon(), Some(0));
}

#[test]
fn verified_artifact_round_trip_is_exact() {
    let directory = tempfile::tempdir().unwrap();
    let artifact_root = directory.path().join("ecological-inputs");
    let state = InitialStateRecipe::Random {
        distribution: DistributionSource::Uniform,
        rng: indexed(903),
    }
    .create(SquareLatticeGeometry::periodic(&[8]).unwrap(), 2)
    .unwrap();
    let persisted = persist_initial_state(&artifact_root, &state).unwrap();
    assert_eq!(persisted.disposition(), ArtifactDisposition::Created);
    assert_eq!(persisted.descriptor().format(), INITIAL_STATE_FORMAT);
    let loaded = load_verified_initial_state(&artifact_root, persisted.descriptor()).unwrap();
    assert_eq!(loaded.space().data(), state.space().data());
    assert_eq!(loaded.counts(), state.counts());
    assert_eq!(loaded.rng_config(), state.rng_config());

    let reused = persist_initial_state(&artifact_root, &state).unwrap();
    assert_eq!(reused.disposition(), ArtifactDisposition::Reused);
    assert_eq!(reused.descriptor(), persisted.descriptor());
    assert_eq!(state.rng_config().unwrap().seed(), 903);
    assert_eq!(
        state.rng_config().unwrap().method(),
        RngMethod::IndexedSplitMix64
    );
    let document: serde_json::Value =
        serde_json::from_slice(&state.to_json_bytes().unwrap()).unwrap();
    assert_eq!(document["format"], "ecological.initial-state.v2");
    assert_eq!(document["rng"]["seed"], 903);
    assert!(document.get("rng_record").is_none());
}

#[test]
fn artifact_corruption_is_rejected_before_initial_state_decoding() {
    let directory = tempfile::tempdir().unwrap();
    let artifact_root = directory.path().join("ecological-inputs");
    let state = InitialStateRecipe::BalancedUniform { rng: indexed(904) }
        .create(SquareLatticeGeometry::periodic(&[8]).unwrap(), 2)
        .unwrap();
    let persisted = persist_initial_state(&artifact_root, &state).unwrap();
    fs::write(
        artifact_root.join(persisted.descriptor().path()),
        b"corrupt",
    )
    .unwrap();

    assert!(matches!(
        load_verified_initial_state(&artifact_root, persisted.descriptor()),
        Err(InitialStateError::ArtifactLoad(
            ArtifactLoadError::DigestMismatch { .. }
        ))
    ));
}

#[test]
fn portable_reference_rejects_root_escape_without_touching_the_target() {
    let directory = tempfile::tempdir().unwrap();
    let artifact_root = directory.path().join("ecological-inputs");
    let state = InitialStateRecipe::BalancedUniform { rng: indexed(905) }
        .create(SquareLatticeGeometry::periodic(&[8]).unwrap(), 2)
        .unwrap();
    let persisted = persist_initial_state(&artifact_root, &state).unwrap();
    let reference =
        InitialStateArtifactReference::new(artifact_root.clone(), persisted.descriptor().clone());
    assert_eq!(reference.artifact_root(), artifact_root);

    let mut value: serde_json::Value =
        serde_json::from_slice(&reference.to_json_bytes().unwrap()).unwrap();
    assert_eq!(
        value["format"],
        "ecological.initial-state-artifact-reference.v2"
    );
    assert!(value.get("artifact_root").is_some());
    assert!(value.get("execution_directory").is_none());
    value["descriptor"]["path"] = "../outside.json".into();
    let escaped =
        InitialStateArtifactReference::from_json_bytes(&serde_json::to_vec(&value).unwrap())
            .unwrap();
    assert!(matches!(
        escaped.resolve(),
        Err(InitialStateError::ArtifactLoad(
            ArtifactLoadError::InvalidDescriptor { .. }
        ))
    ));
}
