use ecological_state_toolkit::artifact::ArtifactDisposition;
use ecological_state_toolkit::interaction::{
    DiagonalPolicy, InteractionArtifactReference, InteractionMatrix, InteractionMatrixRecipe,
    InteractionProvenance, InteractionSourceKind, InteractionTransformation, MatrixNormalization,
    SignStructure, load_verified_interaction_matrix, persist_truncated_svd_series,
};
use physics_in_parallel::prelude::basic::{Backend, Matrix, ResolvedRng, RngMethod};

fn seeded(seed: u64) -> ResolvedRng {
    ResolvedRng::new(seed, RngMethod::IndexedSplitMix64)
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
    let source = InteractionMatrix::from_rows(vec![vec![1.0, 4.0], vec![-2.0, 3.0]]).unwrap();
    assert_eq!(source.species(), 2);
    let mut products = [0.0; 4];
    source
        .mul_vectors_into(&[1.0, 1.0, 2.0, -1.0], &mut products)
        .unwrap();
    assert_eq!(products, [5.0, 1.0, -2.0, -7.0]);
    let antisymmetric = source.antisymmetrize().unwrap();

    assert_eq!(source.coefficient(0, 1), 4.0);
    assert_eq!(antisymmetric.coefficient(0, 0), 0.0);
    assert_eq!(antisymmetric.coefficient(0, 1), 6.0);
    assert_eq!(antisymmetric.coefficient(1, 0), -6.0);

    let scaled = antisymmetric.scale(2.0).unwrap();
    assert_eq!(scaled.coefficient(0, 1), 12.0);

    let absolute = scaled.abs().unwrap();
    assert_eq!(absolute.coefficient(0, 1), 12.0);
    assert_eq!(absolute.coefficient(1, 0), 12.0);
    assert_eq!(scaled.coefficient(1, 0), -12.0);

    let lower_clamped = scaled.clamp_min(0.0).unwrap();
    assert_eq!(lower_clamped.coefficient(0, 0), 0.0);
    assert_eq!(lower_clamped.coefficient(0, 1), 12.0);
    assert_eq!(lower_clamped.coefficient(1, 0), 0.0);
    assert_eq!(scaled.coefficient(1, 0), -12.0);

    let clamped = lower_clamped.clamp_max(10.0).unwrap();
    assert_eq!(clamped.coefficient(0, 0), 0.0);
    assert_eq!(clamped.coefficient(0, 1), 10.0);
    assert_eq!(clamped.coefficient(1, 0), 0.0);
    assert_eq!(lower_clamped.coefficient(0, 1), 12.0);
    assert!(clamped.ensure_max_abs_at_most(10.0).is_ok());
    assert!(matches!(
        lower_clamped.ensure_max_abs_at_most(11.0),
        Err(ecological_state_toolkit::interaction::InteractionMatrixError::MaximumAbsoluteEntryExceeded {
            threshold: 11.0,
            maximum: 12.0,
        })
    ));

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
fn truncated_svd_series_reuses_one_ordered_spectrum_and_preserves_authored_rank_order() {
    let source = InteractionMatrix::from_rows(vec![
        vec![3.0, 0.0, 0.0],
        vec![0.0, 2.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ])
    .unwrap();
    let series = source.truncated_svd_series(&[2, 1, 3]).unwrap();

    assert_eq!(series.singular_values(), &[3.0, 2.0, 1.0]);
    assert_eq!(
        series
            .approximations()
            .iter()
            .map(|approximation| approximation.retained_rank())
            .collect::<Vec<_>>(),
        [2, 1, 3]
    );
    let rank_two = series.approximation(2).unwrap();
    assert_eq!(rank_two.matrix().coefficient(0, 0), 3.0);
    assert_eq!(rank_two.matrix().coefficient(1, 1), 2.0);
    assert_eq!(rank_two.matrix().coefficient(2, 2), 0.0);
    assert!((rank_two.retained_spectral_energy() - 13.0 / 14.0).abs() < 1.0e-12);
    assert!((rank_two.relative_reconstruction_error() - (1.0_f64 / 14.0).sqrt()).abs() < 1.0e-12);
    assert!(matches!(
        rank_two.matrix().provenance(),
        InteractionProvenance::Derived {
            transformation: InteractionTransformation::TruncatedSvd {
                retained_rank: 2,
                ..
            },
            ..
        }
    ));
    let full = series.approximation(3).unwrap();
    assert!(full.relative_reconstruction_error() < 1.0e-12);
}

#[test]
fn truncated_svd_series_rejects_empty_invalid_and_duplicate_ranks() {
    let source = InteractionMatrix::from_rows(vec![vec![1.0, 0.0], vec![0.0, 1.0]]).unwrap();
    assert!(matches!(
        source.truncated_svd_series(&[]),
        Err(ecological_state_toolkit::interaction::InteractionMatrixError::EmptyRetainedRanks)
    ));
    assert!(matches!(
        source.truncated_svd_series(&[0]),
        Err(
            ecological_state_toolkit::interaction::InteractionMatrixError::InvalidRetainedRank {
                rank: 0,
                species: 2
            }
        )
    ));
    assert!(matches!(
        source.truncated_svd_series(&[1, 1]),
        Err(
            ecological_state_toolkit::interaction::InteractionMatrixError::DuplicateRetainedRank {
                rank: 1
            }
        )
    ));
}

#[test]
fn truncated_svd_series_publishes_each_reconstruction_as_verified_json() {
    let root = tempfile::tempdir().unwrap();
    let source = InteractionMatrix::from_rows(vec![
        vec![3.0, 0.0, 0.0],
        vec![0.0, 2.0, 0.0],
        vec![0.0, 0.0, 1.0],
    ])
    .unwrap();
    let series = persist_truncated_svd_series(root.path(), &source, &[1, 2]).unwrap();

    assert_eq!(series.approximations().len(), 2);
    for member in series.approximations() {
        assert_eq!(
            member.persisted().descriptor().path().extension().unwrap(),
            "json"
        );
        let restored =
            load_verified_interaction_matrix(root.path(), member.persisted().descriptor()).unwrap();
        assert_eq!(restored.values(), member.approximation().matrix().values());
    }
}

#[test]
fn matrix_constructors_infer_species_and_reject_non_square_inputs() {
    let matrix = InteractionMatrix::from_matrix(
        Matrix::from_values(2, 2, Backend::Dense, vec![1.0, 2.0, 3.0, 4.0]).unwrap(),
    )
    .unwrap();
    assert_eq!(matrix.species(), 2);

    let error = InteractionMatrix::from_matrix(
        Matrix::from_values(2, 3, Backend::Dense, vec![0.0; 6]).unwrap(),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ecological_state_toolkit::interaction::InteractionMatrixError::NonSquare {
            rows: 2,
            columns: 3
        }
    ));
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
fn zero_connectance_removes_every_structured_off_diagonal_pair() {
    let recipes = [
        InteractionMatrixRecipe::CorrelatedGaussian {
            mean: 0.5,
            standard_deviation: 1.0,
            reciprocal_correlation: -0.25,
            connectance: 0.0,
            diagonal: DiagonalPolicy::Constant(-1.0),
            normalization: MatrixNormalization::SqrtSpecies,
            rng: seeded(92),
        },
        InteractionMatrixRecipe::SignStructuredGaussian {
            structure: SignStructure::Mutualism,
            scale: 1.0,
            connectance: 0.0,
            diagonal: DiagonalPolicy::Constant(-1.0),
            normalization: MatrixNormalization::SqrtSpecies,
            rng: seeded(93),
        },
    ];

    for recipe in recipes {
        let matrix = recipe.generate(6).unwrap();
        for row in 0..6 {
            for column in 0..6 {
                let expected = if row == column { -1.0 } else { 0.0 };
                assert_eq!(matrix.coefficient(row, column), expected);
            }
        }
    }
}

#[test]
fn derived_provenance_survives_verified_artifact_round_trip() {
    let directory = tempfile::tempdir().unwrap();
    let artifact_root = directory.path().join("ecological-inputs");
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
        ecological_state_toolkit::interaction::persist_interaction_matrix(&artifact_root, &matrix)
            .unwrap();
    let loaded = ecological_state_toolkit::interaction::load_verified_interaction_matrix(
        &artifact_root,
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
    assert_eq!(loaded.generator_rng_config(), matrix.generator_rng_config());
    assert_eq!(loaded.generator_rng_config().unwrap().seed(), 104);
    assert_eq!(
        loaded.generator_rng_config().unwrap().method(),
        RngMethod::IndexedSplitMix64
    );

    let reused =
        ecological_state_toolkit::interaction::persist_interaction_matrix(&artifact_root, &matrix)
            .unwrap();
    assert_eq!(reused.disposition(), ArtifactDisposition::Reused);
    assert_eq!(reused.descriptor(), persisted.descriptor());

    let reference =
        InteractionArtifactReference::new(artifact_root.clone(), persisted.into_descriptor());
    assert_eq!(reference.artifact_root(), artifact_root);
    let reference_json: serde_json::Value =
        serde_json::from_slice(&reference.to_json_bytes().unwrap()).unwrap();
    assert_eq!(
        reference_json["format"],
        "ecological.interaction-artifact-reference.v2"
    );
    assert!(reference_json.get("artifact_root").is_some());
    assert!(reference_json.get("execution_directory").is_none());
    let decoded_reference = InteractionArtifactReference::from_json_bytes(
        &serde_json::to_vec(&reference_json).unwrap(),
    )
    .unwrap();
    let reference_loaded = decoded_reference.resolve().unwrap();
    for row in 0..matrix.species() {
        for column in 0..matrix.species() {
            assert_eq!(
                reference_loaded.coefficient(row, column),
                matrix.coefficient(row, column)
            );
        }
    }
}

#[test]
fn obsolete_application_specific_recipe_is_rejected() {
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

    let matrix = InteractionMatrix::from_rows(vec![vec![1.0]]).unwrap();
    assert!(matrix.scale(f64::INFINITY).is_err());
    assert!(matrix.clamp_min(f64::NAN).is_err());
    assert!(matrix.clamp_max(f64::INFINITY).is_err());
    assert!(matrix.ensure_max_abs_at_most(-1.0).is_err());
    assert!(matrix.ensure_max_abs_at_most(f64::INFINITY).is_err());
    assert!(matrix.normalize(-1.0).is_err());
    assert!(matrix.normalize(f64::NAN).is_err());
}
