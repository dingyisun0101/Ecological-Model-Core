use ecological_model_core::interaction::{
    DiagonalPolicy, InteractionMatrixRecipe, MatrixNormalization, SignStructure,
};
use physics_in_parallel::rng::{IndexedRng, RngConfig};
use scientific_workflow::execution::ExecutionScope;

#[test]
fn antisymmetric_recipe_matches_dispatcher_coordinates_exactly() {
    let config = RngConfig::new(Some(17), None, None);
    let matrix = InteractionMatrixRecipe::AntisymmetricGaussian {
        scale: 2.0,
        normalization: MatrixNormalization::SqrtSpecies,
        rng: config,
    }
    .generate(7)
    .unwrap();
    let rng = IndexedRng::new(config).unwrap();
    for row in 0..7 {
        assert_eq!(matrix.coefficient(row, row), 0.0);
        for column in (row + 1)..7 {
            let expected =
                rng.standard_normal(0, 7, row as u64, column as u64) * 2.0 / (7.0_f64).sqrt();
            assert_eq!(matrix.coefficient(row, column), expected);
            assert_eq!(matrix.coefficient(column, row), -expected);
        }
    }
}

#[test]
fn fully_sampled_gaussian_can_be_antisymmetrized_from_the_same_realization() {
    let species = 6;
    let source = InteractionMatrixRecipe::IndependentGaussian {
        mean: 0.0,
        standard_deviation: 1.0,
        connectance: 1.0,
        diagonal: DiagonalPolicy::Sampled,
        normalization: MatrixNormalization::None,
        rng: RngConfig::new(Some(29), None, None),
    }
    .generate(species)
    .unwrap();
    assert!((0..species).any(|index| source.coefficient(index, index) != 0.0));

    let antisymmetric = source.antisymmetrized().unwrap();
    for row in 0..species {
        assert_eq!(antisymmetric.coefficient(row, row), 0.0);
        for column in 0..species {
            assert_eq!(
                antisymmetric.coefficient(row, column),
                (source.coefficient(row, column) - source.coefficient(column, row)) / 2.0
            );
            assert_eq!(
                antisymmetric.coefficient(row, column),
                -antisymmetric.coefficient(column, row)
            );
        }
    }
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
            diagonal: ecological_model_core::interaction::DiagonalPolicy::Zero,
            normalization: MatrixNormalization::SqrtSpecies,
            rng: RngConfig::new(Some(91), None, None),
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
fn generated_provenance_survives_verified_artifact_round_trip() {
    let directory = tempfile::tempdir().unwrap();
    let scope = ExecutionScope::create_named(directory.path(), "execution").unwrap();
    let matrix = InteractionMatrixRecipe::AntisymmetricGaussian {
        scale: 1.0,
        normalization: MatrixNormalization::SqrtSpecies,
        rng: RngConfig::new(Some(104), None, None),
    }
    .generate(4)
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
    assert_eq!(
        loaded.provenance().generator(),
        matrix.provenance().generator()
    );
}
