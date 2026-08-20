use ecological_model_core::interaction::{
    DiagonalPolicy, InteractionMatrix, InteractionMatrixRecipe, InteractionProvenance,
    InteractionSourceKind, InteractionTransformation, MatrixNormalization, SignStructure,
};
use physics_in_parallel::rng::RngConfig;
use scientific_workflow::execution::ExecutionScope;

fn seeded(seed: u64) -> RngConfig {
    RngConfig::new(Some(seed), None)
}

#[test]
fn independent_random_recipes_fill_every_entry_reproducibly() {
    let uniform_recipe = InteractionMatrixRecipe::RandomUniform {
        minimum: -3.0,
        maximum: 2.0,
        rng: seeded(17),
    };
    let first = uniform_recipe.generate(7).unwrap();
    let second = uniform_recipe.generate(7).unwrap();
    for row in 0..7 {
        for column in 0..7 {
            let value = first.coefficient(row, column);
            assert!((-3.0..2.0).contains(&value));
            assert_eq!(value, second.coefficient(row, column));
        }
    }

    let gaussian_recipe = InteractionMatrixRecipe::RandomGaussian {
        mean: 1.5,
        standard_deviation: 0.75,
        rng: seeded(29),
    };
    let first = gaussian_recipe.generate(6).unwrap();
    let second = gaussian_recipe.generate(6).unwrap();
    for row in 0..6 {
        for column in 0..6 {
            assert_eq!(
                first.coefficient(row, column),
                second.coefficient(row, column)
            );
        }
    }
    assert!((0..6).any(|index| first.coefficient(index, index) != 0.0));

    let constant = InteractionMatrixRecipe::RandomGaussian {
        mean: -2.0,
        standard_deviation: 0.0,
        rng: seeded(30),
    }
    .generate(3)
    .unwrap();
    assert_eq!(constant.values().max_abs_real(), 2.0);
    for row in 0..3 {
        for column in 0..3 {
            assert_eq!(constant.coefficient(row, column), -2.0);
        }
    }
}

#[test]
fn transformations_compose_pip_matrix_operations_without_mutating_source() {
    let source = InteractionMatrix::from_rows(vec![vec![1.0, 4.0], vec![-2.0, 3.0]], 2).unwrap();
    let antisymmetric = source.antisymmetrize().unwrap();

    assert_eq!(source.coefficient(0, 1), 4.0);
    assert_eq!(antisymmetric.coefficient(0, 0), 0.0);
    assert_eq!(antisymmetric.coefficient(0, 1), 6.0);
    assert_eq!(antisymmetric.coefficient(1, 0), -6.0);

    let scaled = antisymmetric.scale(2.0).unwrap();
    assert_eq!(scaled.coefficient(0, 1), 12.0);

    let normalized = scaled.normalize(3.0).unwrap();
    assert_eq!(normalized.values().max_abs_real(), 3.0);
    assert_eq!(normalized.coefficient(0, 1), 3.0);
    assert_eq!(normalized.coefficient(1, 0), -3.0);

    let already_bounded = normalized.normalize(10.0).unwrap();
    assert_eq!(already_bounded.coefficient(0, 1), 3.0);
    assert_eq!(
        already_bounded.provenance().kind(),
        InteractionSourceKind::Derived
    );
}

#[test]
fn transformation_provenance_records_the_complete_chain() {
    let source = InteractionMatrixRecipe::RandomGaussian {
        mean: 0.0,
        standard_deviation: 1.0,
        rng: seeded(41),
    }
    .generate(3)
    .unwrap();
    let transformed = source
        .antisymmetrize()
        .unwrap()
        .scale(0.5)
        .unwrap()
        .normalize(1.0)
        .unwrap();

    let InteractionProvenance::Derived {
        source,
        transformation:
            InteractionTransformation::Normalize {
                threshold,
                maximum: _,
                scalar: _,
            },
    } = transformed.provenance()
    else {
        panic!("normalization provenance was not recorded");
    };
    assert_eq!(*threshold, 1.0);
    assert!(matches!(
        source.as_ref(),
        InteractionProvenance::Derived {
            transformation: InteractionTransformation::Scale { scalar: 0.5 },
            ..
        }
    ));
    assert!(transformed.provenance().generator().is_some());
}

#[test]
fn sign_structured_recipes_enforce_ecological_pair_signs() {
    for structure in [
        SignStructure::Competition,
        SignStructure::Mutualism,
        SignStructure::ConsumerResource,
    ] {
        let matrix = InteractionMatrixRecipe::SignStructuredGaussian {
            structure,
            scale: 1.0,
            connectance: 1.0,
            diagonal: DiagonalPolicy::Zero,
            normalization: MatrixNormalization::SqrtSpecies,
            rng: seeded(91),
        }
        .generate(5)
        .unwrap();
        for row in 0..5 {
            for column in (row + 1)..5 {
                let pair = (
                    matrix.coefficient(row, column),
                    matrix.coefficient(column, row),
                );
                match structure {
                    SignStructure::Competition => assert!(pair.0 <= 0.0 && pair.1 <= 0.0),
                    SignStructure::Mutualism => assert!(pair.0 >= 0.0 && pair.1 >= 0.0),
                    SignStructure::ConsumerResource => assert!(pair.0 * pair.1 <= 0.0),
                }
            }
        }
    }
}

#[test]
fn derived_provenance_survives_verified_artifact_round_trip() {
    let directory = tempfile::tempdir().unwrap();
    let scope = ExecutionScope::create_named(directory.path(), "execution").unwrap();
    let matrix = InteractionMatrixRecipe::RandomGaussian {
        mean: 0.0,
        standard_deviation: 1.0,
        rng: seeded(104),
    }
    .generate(4)
    .unwrap()
    .antisymmetrize()
    .unwrap()
    .scale(0.25)
    .unwrap();
    let persisted =
        ecological_model_core::interaction::persist_interaction_matrix(&scope, &matrix).unwrap();
    let loaded = ecological_model_core::interaction::load_verified_interaction_matrix(
        scope.directory(),
        persisted.descriptor(),
    )
    .unwrap();
    for row in 0..4 {
        for column in 0..4 {
            assert_eq!(
                loaded.coefficient(row, column),
                matrix.coefficient(row, column)
            );
        }
    }
    assert_eq!(loaded.provenance(), matrix.provenance());
}

#[test]
fn obsolete_dispatcher_specific_recipe_is_rejected() {
    let legacy = serde_json::json!({
        "family": "antisymmetric_gaussian",
        "scale": 1.0
    });
    assert!(serde_json::from_value::<InteractionMatrixRecipe>(legacy).is_err());
}

#[test]
fn recipes_and_transformations_reject_invalid_numeric_parameters() {
    assert!(
        InteractionMatrixRecipe::RandomUniform {
            minimum: 1.0,
            maximum: 1.0,
            rng: seeded(1),
        }
        .validate(2)
        .is_err()
    );
    assert!(
        InteractionMatrixRecipe::RandomGaussian {
            mean: 0.0,
            standard_deviation: -1.0,
            rng: seeded(1),
        }
        .validate(2)
        .is_err()
    );

    let matrix = InteractionMatrix::from_rows(vec![vec![1.0]], 1).unwrap();
    assert!(matrix.scale(f64::INFINITY).is_err());
    assert!(matrix.normalize(-1.0).is_err());
    assert!(matrix.normalize(f64::NAN).is_err());
}
