use ecological_state_toolkit::terminal_state::{
    StopReason, TerminalClassification, TerminationSignal,
};
use ecological_state_toolkit::trajectory::{
    AbundanceView, DetectionPolicy, EquilibriumEvidence, EquilibriumPolicy, PeriodicOrbitPolicy,
    ResidualTolerance, TerminalPolicy, TrajectoryObservation, TrajectoryObservationPolicy,
    TrajectoryObserver,
};

fn observation<'a>(iteration: u64, abundance: &'a [f64]) -> TrajectoryObservation<'a> {
    TrajectoryObservation {
        iteration,
        physical_time: Some(iteration as f64),
        abundance: AbundanceView::Continuous(abundance),
        detector_observable: None,
        equilibrium_evidence: EquilibriumEvidence::Unavailable,
    }
}

#[test]
fn disabled_constructs_no_observer() {
    assert!(
        TrajectoryObserver::from_policy(TrajectoryObservationPolicy::Disabled)
            .unwrap()
            .is_none()
    );
}

#[test]
fn terminal_only_forces_an_off_cadence_final_sample() {
    let mut observer = TrajectoryObserver::from_policy(TrajectoryObservationPolicy::TerminalOnly(
        TerminalPolicy {
            sample_interval_iterations: 2,
            trailing_window_samples: 2,
        },
    ))
    .unwrap()
    .unwrap();
    let states = [[3.0, 1.0], [2.0, 2.0], [1.0, 3.0], [4.0, 1.0]];
    for (iteration, state) in states.iter().enumerate() {
        assert!(
            observer
                .observe(observation(iteration as u64, state))
                .unwrap()
                .is_none()
        );
    }
    let terminal = observer
        .finish(observation(3, &states[3]), StopReason::MaximumIterations)
        .unwrap();
    assert_eq!(
        terminal.classification(),
        TerminalClassification::TrailingAverage
    );
    assert_eq!(terminal.first_sample_iteration(), 2);
    assert_eq!(terminal.last_sample_iteration(), 3);
    assert_eq!(terminal.composition(), &[0.525, 0.475]);
}

#[test]
fn staged_equilibrium_requires_authoritative_evidence() {
    let policy = DetectionPolicy {
        terminal: TerminalPolicy {
            sample_interval_iterations: 1,
            trailing_window_samples: 3,
        },
        start_after_iteration: 0,
        equilibrium: Some(EquilibriumPolicy {
            base_window_samples: 2,
            confirmation_window_multipliers: vec![1, 2],
            maximum_observable_distance: 1.0e-12,
            maximum_relative_mass_range: Some(0.0),
            support_threshold: 0.0,
            residual_tolerance: ResidualTolerance {
                absolute: 1.0e-8,
                relative: 1.0e-6,
            },
        }),
        periodic_orbit: None,
        detect_absorbing_state: true,
    };
    let mut observer = TrajectoryObserver::from_policy(TrajectoryObservationPolicy::Detect(policy))
        .unwrap()
        .unwrap();
    let state = [1.0, 1.0];
    let mut signal = None;
    for iteration in 0..6 {
        let mut value = observation(iteration, &state);
        value.equilibrium_evidence = EquilibriumEvidence::MaximumScaledResidual { value: 0.5 };
        signal = observer.observe(value).unwrap();
    }
    assert!(matches!(signal, Some(TerminationSignal::Equilibrium(_))));
}

#[test]
fn recurrent_nonconstant_states_form_a_periodic_orbit() {
    let policy = DetectionPolicy {
        terminal: TerminalPolicy {
            sample_interval_iterations: 1,
            trailing_window_samples: 4,
        },
        start_after_iteration: 0,
        equilibrium: None,
        periodic_orbit: Some(PeriodicOrbitPolicy {
            minimum_period_samples: 2,
            maximum_period_samples: 2,
            repeated_cycles: 2,
            maximum_recurrence_distance: 0.0,
            minimum_orbit_amplitude: 0.1,
        }),
        detect_absorbing_state: false,
    };
    let mut observer = TrajectoryObserver::from_policy(TrajectoryObservationPolicy::Detect(policy))
        .unwrap()
        .unwrap();
    let states = [[0.9, 0.1], [0.1, 0.9]];
    let mut signal = None;
    for iteration in 0..5 {
        signal = observer
            .observe(observation(iteration, &states[iteration as usize % 2]))
            .unwrap();
    }
    assert!(matches!(signal, Some(TerminationSignal::PeriodicOrbit(_))));
}

#[test]
fn explicit_absorption_yields_the_exact_count_composition() {
    let policy = DetectionPolicy {
        terminal: TerminalPolicy {
            sample_interval_iterations: 1,
            trailing_window_samples: 2,
        },
        start_after_iteration: 0,
        equilibrium: None,
        periodic_orbit: None,
        detect_absorbing_state: true,
    };
    let mut observer = TrajectoryObserver::from_policy(TrajectoryObservationPolicy::Detect(policy))
        .unwrap()
        .unwrap();
    let counts = [1_usize, 3];
    let value = TrajectoryObservation {
        iteration: 4,
        physical_time: None,
        abundance: AbundanceView::Counts(&counts),
        detector_observable: None,
        equilibrium_evidence: EquilibriumEvidence::AbsorbingState,
    };
    let signal = observer.observe(value).unwrap().unwrap();
    let terminal = observer
        .finish(value, StopReason::Detected(signal))
        .unwrap();
    assert_eq!(
        terminal.classification(),
        TerminalClassification::AbsorbingState
    );
    assert_eq!(terminal.composition(), &[0.25, 0.75]);
    let decoded = ecological_state_toolkit::terminal_state::TerminalState::from_json_bytes(
        &terminal.to_json_bytes().unwrap(),
    )
    .unwrap();
    assert_eq!(decoded, terminal);
}
