use ecological_model_core::initial_state::{
    DistributionSource, INITIAL_STATE_FORMAT, InitialStateRecipe, InitializationMethod,
    load_verified_initial_state, persist_initial_state,
};
use physics_in_parallel::rng::RngConfig;
use physics_in_parallel::space::discrete::square_lattice::SquareLatticeConfig;
use scientific_workflow::artifact::ArtifactDisposition;
use scientific_workflow::execution::ExecutionScope;

#[test]
fn recipe_is_reproducible_and_counts_are_exact() {
    let lattice = SquareLatticeConfig::periodic(&[9, 9]);
    let recipe = InitialStateRecipe::CenteredSeed {
        distribution: DistributionSource::Inline {
            weights: vec![0.7, 0.3],
        },
        seed_taxon: 1,
        seed_radius: 1,
        rng: RngConfig::new(Some(42), None, None),
    };
    let first = recipe.clone().create(lattice.clone(), 2, None).unwrap();
    let second = recipe.create(lattice, 2, None).unwrap();
    assert_eq!(first.space().data(), second.space().data());
    assert_eq!(first.counts().iter().sum::<usize>(), 81);
    assert_eq!(first.method(), InitializationMethod::CenteredSeed);
    assert_eq!(first.seed_taxon(), Some(1));
}

#[test]
fn dominant_recipe_uses_first_maximal_taxon() {
    let state = InitialStateRecipe::CenteredDominantSeed {
        distribution: DistributionSource::Inline {
            weights: vec![0.4, 0.4, 0.2],
        },
        seed_radius: 1,
        rng: RngConfig::new(Some(7), None, None),
    }
    .create(SquareLatticeConfig::periodic(&[9]), 3, None)
    .unwrap();
    assert_eq!(state.seed_taxon(), Some(0));
    assert_eq!(state.counts()[0], 3);
}

#[test]
fn verified_artifact_round_trip_is_exact() {
    let directory = tempfile::tempdir().unwrap();
    let scope = ExecutionScope::create_named(directory.path(), "execution").unwrap();
    let state = InitialStateRecipe::Random {
        distribution: DistributionSource::Uniform,
        rng: RngConfig::new(Some(903), None, None),
    }
    .create(SquareLatticeConfig::periodic(&[8]), 2, None)
    .unwrap();
    let persisted = persist_initial_state(&scope, &state).unwrap();
    assert_eq!(persisted.disposition(), ArtifactDisposition::Created);
    assert_eq!(persisted.descriptor().format(), INITIAL_STATE_FORMAT);
    let loaded = load_verified_initial_state(scope.directory(), persisted.descriptor()).unwrap();
    assert_eq!(loaded.space().data(), state.space().data());
    assert_eq!(loaded.counts(), state.counts());
}
