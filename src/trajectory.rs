//! Allocation-conscious observation of ecological trajectories.

use crate::terminal_state::{
    AbsorptionDiagnostics, EquilibriumDiagnostics, PeriodicOrbitDiagnostics, StopReason,
    TerminalClassification, TerminalState, TerminalStateError, TerminationSignal,
};
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ResidualTolerance {
    pub absolute: f64,
    pub relative: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct TerminalPolicy {
    pub sample_interval_iterations: u64,
    pub trailing_window_samples: usize,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct EquilibriumPolicy {
    pub base_window_samples: usize,
    pub confirmation_window_multipliers: Vec<usize>,
    pub maximum_observable_distance: f64,
    pub maximum_relative_mass_range: Option<f64>,
    pub support_threshold: f64,
    pub residual_tolerance: ResidualTolerance,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct PeriodicOrbitPolicy {
    pub minimum_period_samples: usize,
    pub maximum_period_samples: usize,
    pub repeated_cycles: usize,
    pub maximum_recurrence_distance: f64,
    pub minimum_orbit_amplitude: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct DetectionPolicy {
    pub terminal: TerminalPolicy,
    #[serde(default)]
    pub start_after_iteration: u64,
    pub equilibrium: Option<EquilibriumPolicy>,
    pub periodic_orbit: Option<PeriodicOrbitPolicy>,
    #[serde(default = "default_true")]
    pub detect_absorbing_state: bool,
}

const fn default_true() -> bool {
    true
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "mode", content = "policy", rename_all = "snake_case")]
pub enum TrajectoryObservationPolicy {
    Disabled,
    TerminalOnly(TerminalPolicy),
    Detect(DetectionPolicy),
}

#[derive(Clone, Copy, Debug)]
pub enum AbundanceView<'a> {
    Continuous(&'a [f64]),
    Counts(&'a [usize]),
}

#[derive(Clone, Copy, Debug)]
pub enum EquilibriumEvidence<'a> {
    Unavailable,
    Residual { values: &'a [f64] },
    MaximumScaledResidual { value: f64 },
    AbsorbingState,
}

#[derive(Clone, Copy, Debug)]
pub struct TrajectoryObservation<'a> {
    pub iteration: u64,
    pub physical_time: Option<f64>,
    pub abundance: AbundanceView<'a>,
    pub detector_observable: Option<&'a [f64]>,
    pub equilibrium_evidence: EquilibriumEvidence<'a>,
}

#[derive(Default)]
struct Sample {
    iteration: u64,
    physical_time: Option<f64>,
    composition: Vec<f64>,
    observable: Option<Vec<f64>>,
    mass: f64,
    scaled_residual: Option<f64>,
    absorbing: bool,
}

impl Sample {
    fn detector_observable(&self) -> &[f64] {
        self.observable.as_deref().unwrap_or(&self.composition)
    }
}

struct SampleArena {
    slots: Vec<Sample>,
    head: usize,
    len: usize,
}

impl SampleArena {
    fn new(capacity: usize) -> Self {
        Self {
            slots: (0..capacity).map(|_| Sample::default()).collect(),
            head: 0,
            len: 0,
        }
    }

    fn push_slot(&mut self) -> &mut Sample {
        let index = if self.len < self.slots.len() {
            let index = (self.head + self.len) % self.slots.len();
            self.len += 1;
            index
        } else {
            let index = self.head;
            self.head = (self.head + 1) % self.slots.len();
            index
        };
        &mut self.slots[index]
    }

    fn get(&self, offset: usize) -> &Sample {
        &self.slots[(self.head + offset) % self.slots.len()]
    }

    fn last(&self) -> Option<&Sample> {
        self.len.checked_sub(1).map(|index| self.get(index))
    }

    fn take_last_composition(&mut self) -> Option<Vec<f64>> {
        let index = self.len.checked_sub(1)?;
        let slot = (self.head + index) % self.slots.len();
        Some(std::mem::take(&mut self.slots[slot].composition))
    }
}

/// A bounded observer. Construct it with [`TrajectoryObserver::from_policy`]
/// so disabled observation allocates nothing.
pub struct TrajectoryObserver {
    policy: ActivePolicy,
    arena: SampleArena,
    abundance_dimension: Option<usize>,
    observable_dimension: Option<usize>,
    uses_distinct_observable: Option<bool>,
    last_seen_iteration: Option<u64>,
    last_seen_time: Option<f64>,
    equilibrium_stage: usize,
    equilibrium_stage_samples: usize,
    scratch: Vec<f64>,
}

enum ActivePolicy {
    TerminalOnly(TerminalPolicy),
    Detect(DetectionPolicy),
}

impl TrajectoryObserver {
    /// Returns `None` for `Disabled`, ensuring that mode has no observer or allocation.
    pub fn from_policy(
        policy: TrajectoryObservationPolicy,
    ) -> Result<Option<Self>, TrajectoryObserverError> {
        let active = match policy {
            TrajectoryObservationPolicy::Disabled => return Ok(None),
            TrajectoryObservationPolicy::TerminalOnly(policy) => {
                validate_terminal_policy(&policy)?;
                ActivePolicy::TerminalOnly(policy)
            }
            TrajectoryObservationPolicy::Detect(policy) => {
                validate_detection_policy(&policy)?;
                ActivePolicy::Detect(policy)
            }
        };
        let capacity = required_capacity(&active)?;
        Ok(Some(Self {
            policy: active,
            arena: SampleArena::new(capacity),
            abundance_dimension: None,
            observable_dimension: None,
            uses_distinct_observable: None,
            last_seen_iteration: None,
            last_seen_time: None,
            equilibrium_stage: 0,
            equilibrium_stage_samples: 0,
            scratch: Vec::new(),
        }))
    }

    pub const fn sample_interval_iterations(&self) -> u64 {
        match &self.policy {
            ActivePolicy::TerminalOnly(policy) => policy.sample_interval_iterations,
            ActivePolicy::Detect(policy) => policy.terminal.sample_interval_iterations,
        }
    }

    /// Lets a model avoid calculating a residual on observations that will not be retained.
    pub fn requires_equilibrium_evidence(&self, iteration: u64) -> bool {
        self.is_due(iteration)
            && matches!(&self.policy, ActivePolicy::Detect(policy)
                if iteration >= policy.start_after_iteration && policy.equilibrium.is_some())
    }

    pub fn observe(
        &mut self,
        observation: TrajectoryObservation<'_>,
    ) -> Result<Option<TerminationSignal>, TrajectoryObserverError> {
        self.validate_order(&observation)?;
        let due = self.is_due(observation.iteration);
        self.last_seen_iteration = Some(observation.iteration);
        self.last_seen_time = observation.physical_time.or(self.last_seen_time);
        if !due {
            return Ok(None);
        }
        self.retain(observation)?;
        self.detect()
    }

    /// Completes the observer and always incorporates `final_observation`, even off cadence.
    pub fn finish(
        mut self,
        final_observation: TrajectoryObservation<'_>,
        stop_reason: StopReason,
    ) -> Result<TerminalState, TrajectoryObserverError> {
        let already_retained = self
            .arena
            .last()
            .is_some_and(|sample| sample.iteration == final_observation.iteration);
        if !already_retained {
            self.validate_final_order(&final_observation)?;
            self.retain(final_observation)?;
        }
        let final_sample = self
            .arena
            .last()
            .ok_or(TrajectoryObserverError::NoSamples)?;
        if final_sample.iteration
            != stop_reason
                .signal()
                .map_or(final_sample.iteration, TerminationSignal::iteration)
        {
            return Err(TrajectoryObserverError::StopReasonIterationMismatch);
        }
        let final_iteration = final_sample.iteration;
        let final_time = final_sample.physical_time;

        match &stop_reason {
            StopReason::Detected(TerminationSignal::Equilibrium(_)) => TerminalState::new(
                TerminalClassification::Equilibrium,
                stop_reason,
                final_iteration,
                final_time,
                self.arena
                    .take_last_composition()
                    .expect("final sample was validated"),
                1,
                final_iteration,
                final_iteration,
            ),
            StopReason::Detected(TerminationSignal::AbsorbingState(_)) => TerminalState::new(
                TerminalClassification::AbsorbingState,
                stop_reason,
                final_iteration,
                final_time,
                self.arena
                    .take_last_composition()
                    .expect("final sample was validated"),
                1,
                final_iteration,
                final_iteration,
            ),
            StopReason::Detected(TerminationSignal::PeriodicOrbit(diagnostics)) => {
                let (composition, count, first, last) = self.average_range(
                    diagnostics.first_cycle_iteration,
                    diagnostics.last_cycle_iteration,
                )?;
                TerminalState::new(
                    TerminalClassification::PeriodicOrbit,
                    stop_reason,
                    final_iteration,
                    final_time,
                    composition,
                    count,
                    first,
                    last,
                )
            }
            _ => {
                let count = self.trailing_window_samples().min(self.arena.len);
                let start = self.arena.len - count;
                let first = self.arena.get(start).iteration;
                let last = self.arena.get(self.arena.len - 1).iteration;
                let composition = average_samples(&self.arena, start, self.arena.len);
                TerminalState::new(
                    TerminalClassification::TrailingAverage,
                    stop_reason,
                    final_iteration,
                    final_time,
                    composition,
                    count,
                    first,
                    last,
                )
            }
        }
        .map_err(Into::into)
    }

    fn validate_order(
        &self,
        observation: &TrajectoryObservation<'_>,
    ) -> Result<(), TrajectoryObserverError> {
        if self
            .last_seen_iteration
            .is_some_and(|last| observation.iteration <= last)
        {
            return Err(TrajectoryObserverError::NonIncreasingIteration);
        }
        if observation
            .physical_time
            .is_some_and(|time| !time.is_finite())
            || matches!((self.last_seen_time, observation.physical_time), (Some(last), Some(time)) if time <= last)
        {
            return Err(TrajectoryObserverError::InvalidPhysicalTime);
        }
        Ok(())
    }

    fn validate_final_order(
        &self,
        observation: &TrajectoryObservation<'_>,
    ) -> Result<(), TrajectoryObserverError> {
        if self
            .last_seen_iteration
            .is_some_and(|last| observation.iteration < last)
        {
            return Err(TrajectoryObserverError::NonIncreasingIteration);
        }
        if observation
            .physical_time
            .is_some_and(|time| !time.is_finite())
            || matches!((self.last_seen_iteration, self.last_seen_time, observation.physical_time),
                (Some(last_iteration), Some(last_time), Some(time))
                    if observation.iteration == last_iteration && time != last_time
                        || observation.iteration > last_iteration && time <= last_time)
        {
            return Err(TrajectoryObserverError::InvalidPhysicalTime);
        }
        Ok(())
    }

    fn is_due(&self, iteration: u64) -> bool {
        iteration.is_multiple_of(self.sample_interval_iterations())
    }

    fn trailing_window_samples(&self) -> usize {
        match &self.policy {
            ActivePolicy::TerminalOnly(policy) => policy.trailing_window_samples,
            ActivePolicy::Detect(policy) => policy.terminal.trailing_window_samples,
        }
    }

    fn retain(
        &mut self,
        observation: TrajectoryObservation<'_>,
    ) -> Result<(), TrajectoryObserverError> {
        let abundance_len = abundance_len(observation.abundance);
        require_stable_dimension(&mut self.abundance_dimension, abundance_len, "abundance")?;
        let observable_len = observation
            .detector_observable
            .map_or(abundance_len, <[f64]>::len);
        require_stable_dimension(
            &mut self.observable_dimension,
            observable_len,
            "detector observable",
        )?;
        let uses_distinct = observation.detector_observable.is_some();
        match self.uses_distinct_observable {
            Some(expected) if expected != uses_distinct => {
                return Err(TrajectoryObserverError::DetectorObservableModeChanged);
            }
            None => self.uses_distinct_observable = Some(uses_distinct),
            _ => {}
        }

        let equilibrium = match &self.policy {
            ActivePolicy::Detect(policy)
                if observation.iteration >= policy.start_after_iteration =>
            {
                policy.equilibrium.as_ref()
            }
            _ => None,
        };
        let (mass, scaled_residual, absorbing) = validate_evidence(
            observation.abundance,
            observation.detector_observable,
            observation.equilibrium_evidence,
            equilibrium,
        )?;
        let slot = self.arena.push_slot();
        slot.iteration = observation.iteration;
        slot.physical_time = observation.physical_time;
        slot.mass = mass;
        slot.scaled_residual = scaled_residual;
        slot.absorbing = absorbing;
        normalize_abundance(observation.abundance, &mut slot.composition)?;
        if let Some(values) = observation.detector_observable {
            normalize_continuous(values, slot.observable.get_or_insert_with(Vec::new))?;
        } else {
            slot.observable = None;
        }
        Ok(())
    }

    fn detect(&mut self) -> Result<Option<TerminationSignal>, TrajectoryObserverError> {
        let ActivePolicy::Detect(policy) = &self.policy else {
            return Ok(None);
        };
        let Some(current) = self.arena.last() else {
            return Ok(None);
        };
        if current.iteration < policy.start_after_iteration {
            return Ok(None);
        }
        if policy.detect_absorbing_state && current.absorbing {
            let threshold = policy
                .equilibrium
                .as_ref()
                .map_or(0.0, |value| value.support_threshold);
            let supported_taxa = current
                .composition
                .iter()
                .filter(|value| **value > threshold)
                .count();
            return Ok(Some(TerminationSignal::AbsorbingState(
                AbsorptionDiagnostics {
                    iteration: current.iteration,
                    supported_taxa,
                },
            )));
        }
        if policy.equilibrium.is_some()
            && let Some(signal) = self.detect_equilibrium()?
        {
            return Ok(Some(signal));
        }
        Ok(self.detect_periodic_orbit())
    }

    fn detect_equilibrium(&mut self) -> Result<Option<TerminationSignal>, TrajectoryObserverError> {
        let ActivePolicy::Detect(policy) = &self.policy else {
            return Ok(None);
        };
        let equilibrium = policy
            .equilibrium
            .as_ref()
            .expect("checked by caller")
            .clone();
        if self.arena.len >= 2
            && !same_support(
                self.arena.get(self.arena.len - 2).detector_observable(),
                self.arena.get(self.arena.len - 1).detector_observable(),
                equilibrium.support_threshold,
            )
        {
            self.equilibrium_stage = 0;
            self.equilibrium_stage_samples = 0;
        }
        let required = equilibrium.base_window_samples
            * equilibrium.confirmation_window_multipliers[self.equilibrium_stage];
        self.equilibrium_stage_samples += 1;
        if self.equilibrium_stage_samples < required || self.arena.len < required {
            return Ok(None);
        }
        self.equilibrium_stage_samples = 0;
        let start = self.arena.len - required;
        let maximum_distance = self.maximum_arena_distance_from_mean(start);
        let passed = window_support_is_invariant(&self.arena, start, equilibrium.support_threshold)
            && relative_mass_range(&self.arena, start)
                <= equilibrium
                    .maximum_relative_mass_range
                    .unwrap_or(f64::INFINITY)
            && maximum_distance <= equilibrium.maximum_observable_distance
            && (start..self.arena.len).all(|index| {
                self.arena
                    .get(index)
                    .scaled_residual
                    .is_some_and(|value| value <= 1.0)
            });
        if !passed {
            self.equilibrium_stage = 0;
            return Ok(None);
        }
        self.equilibrium_stage += 1;
        if self.equilibrium_stage < equilibrium.confirmation_window_multipliers.len() {
            return Ok(None);
        }
        let current = self.arena.last().expect("nonempty window");
        let max_residual = (start..self.arena.len)
            .filter_map(|index| self.arena.get(index).scaled_residual)
            .fold(0.0, f64::max);
        Ok(Some(TerminationSignal::Equilibrium(
            EquilibriumDiagnostics {
                iteration: current.iteration,
                completed_windows: self.equilibrium_stage,
                final_window_samples: required,
                maximum_observable_distance: maximum_distance,
                relative_mass_range: relative_mass_range(&self.arena, start),
                maximum_scaled_residual: max_residual,
            },
        )))
    }

    fn maximum_arena_distance_from_mean(&mut self, start: usize) -> f64 {
        fill_mean_observable(&self.arena, start, &mut self.scratch);
        (start..self.arena.len)
            .map(|index| jensen_shannon(self.arena.get(index).detector_observable(), &self.scratch))
            .fold(0.0, f64::max)
    }

    fn detect_periodic_orbit(&self) -> Option<TerminationSignal> {
        let ActivePolicy::Detect(policy) = &self.policy else {
            return None;
        };
        let periodic = policy.periodic_orbit.as_ref()?;
        for period in periodic.minimum_period_samples..=periodic.maximum_period_samples {
            let required = period * periodic.repeated_cycles + 1;
            if self.arena.len < required {
                continue;
            }
            let start = self.arena.len - required;
            let mut maximum_recurrence: f64 = 0.0;
            let mut amplitude: f64 = 0.0;
            for index in (start + period)..self.arena.len {
                maximum_recurrence = maximum_recurrence.max(jensen_shannon(
                    self.arena.get(index).detector_observable(),
                    self.arena.get(index - period).detector_observable(),
                ));
            }
            for index in (start + 1)..(start + period + 1) {
                amplitude = amplitude.max(jensen_shannon(
                    self.arena.get(start).detector_observable(),
                    self.arena.get(index).detector_observable(),
                ));
            }
            if maximum_recurrence <= periodic.maximum_recurrence_distance
                && amplitude >= periodic.minimum_orbit_amplitude
            {
                let current = self.arena.last().expect("required samples");
                return Some(TerminationSignal::PeriodicOrbit(PeriodicOrbitDiagnostics {
                    iteration: current.iteration,
                    period_samples: period,
                    repeated_cycles: periodic.repeated_cycles,
                    first_cycle_iteration: self.arena.get(start).iteration,
                    last_cycle_iteration: current.iteration,
                    maximum_recurrence_distance: maximum_recurrence,
                    orbit_amplitude: amplitude,
                }));
            }
        }
        None
    }

    fn average_range(
        &self,
        first_iteration: u64,
        last_iteration: u64,
    ) -> Result<(Vec<f64>, usize, u64, u64), TrajectoryObserverError> {
        let start = (0..self.arena.len)
            .find(|index| self.arena.get(*index).iteration == first_iteration)
            .ok_or(TrajectoryObserverError::DetectionHistoryUnavailable)?;
        let end = (start..self.arena.len)
            .find(|index| self.arena.get(*index).iteration == last_iteration)
            .ok_or(TrajectoryObserverError::DetectionHistoryUnavailable)?
            + 1;
        Ok((
            average_samples(&self.arena, start, end),
            end - start,
            first_iteration,
            last_iteration,
        ))
    }
}

fn validate_terminal_policy(policy: &TerminalPolicy) -> Result<(), TrajectoryObserverError> {
    if policy.sample_interval_iterations == 0 || policy.trailing_window_samples == 0 {
        return Err(TrajectoryObserverError::InvalidPolicy(
            "terminal sampling values must be positive",
        ));
    }
    Ok(())
}

fn validate_detection_policy(policy: &DetectionPolicy) -> Result<(), TrajectoryObserverError> {
    validate_terminal_policy(&policy.terminal)?;
    if let Some(value) = &policy.equilibrium
        && (value.base_window_samples == 0
            || value.confirmation_window_multipliers.is_empty()
            || value.confirmation_window_multipliers.contains(&0)
            || !valid_nonnegative(value.maximum_observable_distance)
            || !valid_nonnegative(value.support_threshold)
            || value
                .maximum_relative_mass_range
                .is_some_and(|limit| !valid_nonnegative(limit))
            || !valid_positive(value.residual_tolerance.absolute)
            || !valid_nonnegative(value.residual_tolerance.relative))
    {
        return Err(TrajectoryObserverError::InvalidPolicy(
            "invalid equilibrium policy",
        ));
    }
    if let Some(value) = &policy.periodic_orbit
        && (value.minimum_period_samples == 0
            || value.maximum_period_samples < value.minimum_period_samples
            || value.repeated_cycles < 2
            || !valid_nonnegative(value.maximum_recurrence_distance)
            || !valid_positive(value.minimum_orbit_amplitude))
    {
        return Err(TrajectoryObserverError::InvalidPolicy(
            "invalid periodic-orbit policy",
        ));
    }
    Ok(())
}

fn required_capacity(policy: &ActivePolicy) -> Result<usize, TrajectoryObserverError> {
    let terminal = match policy {
        ActivePolicy::TerminalOnly(value) => return Ok(value.trailing_window_samples),
        ActivePolicy::Detect(value) => &value.terminal,
    };
    let ActivePolicy::Detect(detection) = policy else {
        unreachable!()
    };
    let equilibrium = detection
        .equilibrium
        .as_ref()
        .map_or(Some(0), |value| {
            value.base_window_samples.checked_mul(
                value
                    .confirmation_window_multipliers
                    .iter()
                    .copied()
                    .max()
                    .unwrap_or(1),
            )
        })
        .ok_or(TrajectoryObserverError::InvalidPolicy(
            "history capacity overflows usize",
        ))?;
    let periodic = detection
        .periodic_orbit
        .as_ref()
        .map_or(Some(0), |value| {
            value
                .maximum_period_samples
                .checked_mul(value.repeated_cycles)?
                .checked_add(1)
        })
        .ok_or(TrajectoryObserverError::InvalidPolicy(
            "history capacity overflows usize",
        ))?;
    Ok(terminal
        .trailing_window_samples
        .max(equilibrium)
        .max(periodic))
}

fn validate_evidence(
    abundance: AbundanceView<'_>,
    observable: Option<&[f64]>,
    evidence: EquilibriumEvidence<'_>,
    policy: Option<&EquilibriumPolicy>,
) -> Result<(f64, Option<f64>, bool), TrajectoryObserverError> {
    let mass = validate_abundance(abundance)?;
    if let Some(values) = observable {
        validate_continuous(values)?;
    }
    let Some(policy) = policy else {
        return match evidence {
            EquilibriumEvidence::MaximumScaledResidual { value } if !valid_nonnegative(value) => {
                Err(TrajectoryObserverError::InvalidEvidence)
            }
            EquilibriumEvidence::Residual { values }
                if values.iter().any(|value| !value.is_finite()) =>
            {
                Err(TrajectoryObserverError::InvalidEvidence)
            }
            EquilibriumEvidence::AbsorbingState => Ok((mass, None, true)),
            _ => Ok((mass, None, false)),
        };
    };
    match evidence {
        EquilibriumEvidence::Unavailable => Ok((mass, None, false)),
        EquilibriumEvidence::AbsorbingState => Ok((mass, None, true)),
        EquilibriumEvidence::MaximumScaledResidual { value } if valid_nonnegative(value) => {
            Ok((mass, Some(value), false))
        }
        EquilibriumEvidence::MaximumScaledResidual { .. } => {
            Err(TrajectoryObserverError::InvalidEvidence)
        }
        EquilibriumEvidence::Residual { values } => {
            if values.len() != observable.map_or_else(|| abundance_len(abundance), <[f64]>::len)
                || values.iter().any(|value| !value.is_finite())
            {
                return Err(TrajectoryObserverError::InvalidEvidence);
            }
            let scaled = |residual: &f64, state: f64| {
                residual.abs()
                    / (policy.residual_tolerance.absolute
                        + policy.residual_tolerance.relative * state.abs())
            };
            let maximum = if let Some(observed) = observable {
                values
                    .iter()
                    .zip(observed)
                    .map(|(residual, state)| scaled(residual, *state))
                    .fold(0.0, f64::max)
            } else {
                match abundance {
                    AbundanceView::Continuous(observed) => values
                        .iter()
                        .zip(observed)
                        .map(|(residual, state)| scaled(residual, *state))
                        .fold(0.0, f64::max),
                    AbundanceView::Counts(observed) => values
                        .iter()
                        .zip(observed)
                        .map(|(residual, state)| scaled(residual, *state as f64))
                        .fold(0.0, f64::max),
                }
            };
            Ok((mass, Some(maximum), false))
        }
    }
}

fn abundance_len(value: AbundanceView<'_>) -> usize {
    match value {
        AbundanceView::Continuous(v) => v.len(),
        AbundanceView::Counts(v) => v.len(),
    }
}

fn validate_abundance(value: AbundanceView<'_>) -> Result<f64, TrajectoryObserverError> {
    match value {
        AbundanceView::Continuous(values) => validate_continuous(values),
        AbundanceView::Counts(values) => {
            if values.is_empty() {
                return Err(TrajectoryObserverError::InvalidAbundance);
            }
            let total = values
                .iter()
                .try_fold(0usize, |sum, value| sum.checked_add(*value))
                .ok_or(TrajectoryObserverError::InvalidAbundance)?;
            if total == 0 {
                Err(TrajectoryObserverError::InvalidAbundance)
            } else {
                Ok(total as f64)
            }
        }
    }
}

fn validate_continuous(values: &[f64]) -> Result<f64, TrajectoryObserverError> {
    if values.is_empty()
        || values
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(TrajectoryObserverError::InvalidAbundance);
    }
    let total = values.iter().sum::<f64>();
    if valid_positive(total) {
        Ok(total)
    } else {
        Err(TrajectoryObserverError::InvalidAbundance)
    }
}

fn normalize_abundance(
    value: AbundanceView<'_>,
    output: &mut Vec<f64>,
) -> Result<(), TrajectoryObserverError> {
    let total = validate_abundance(value)?;
    output.clear();
    match value {
        AbundanceView::Continuous(values) => {
            output.extend(values.iter().map(|value| value / total))
        }
        AbundanceView::Counts(values) => {
            output.extend(values.iter().map(|value| *value as f64 / total))
        }
    }
    Ok(())
}

fn normalize_continuous(
    values: &[f64],
    output: &mut Vec<f64>,
) -> Result<(), TrajectoryObserverError> {
    let total = validate_continuous(values)?;
    output.clear();
    output.extend(values.iter().map(|value| value / total));
    Ok(())
}

fn require_stable_dimension(
    target: &mut Option<usize>,
    actual: usize,
    name: &'static str,
) -> Result<(), TrajectoryObserverError> {
    match *target {
        Some(expected) if expected != actual => Err(TrajectoryObserverError::DimensionChanged {
            name,
            expected,
            actual,
        }),
        None => {
            *target = Some(actual);
            Ok(())
        }
        _ => Ok(()),
    }
}

fn window_support_is_invariant(arena: &SampleArena, start: usize, threshold: f64) -> bool {
    let reference = arena.get(start).detector_observable();
    (start + 1..arena.len)
        .all(|index| same_support(reference, arena.get(index).detector_observable(), threshold))
}

fn same_support(left: &[f64], right: &[f64], threshold: f64) -> bool {
    left.iter()
        .zip(right)
        .all(|(left, right)| (*left > threshold) == (*right > threshold))
}

fn relative_mass_range(arena: &SampleArena, start: usize) -> f64 {
    let (minimum, maximum) = (start..arena.len).map(|index| arena.get(index).mass).fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(minimum, maximum), value| (minimum.min(value), maximum.max(value)),
    );
    (maximum - minimum) / maximum.max(f64::MIN_POSITIVE)
}

fn fill_mean_observable(arena: &SampleArena, start: usize, output: &mut Vec<f64>) {
    output.clear();
    output.resize(arena.get(start).detector_observable().len(), 0.0);
    for index in start..arena.len {
        for (sum, value) in output
            .iter_mut()
            .zip(arena.get(index).detector_observable())
        {
            *sum += value;
        }
    }
    let count = (arena.len - start) as f64;
    for value in output {
        *value /= count;
    }
}

fn average_samples(arena: &SampleArena, start: usize, end: usize) -> Vec<f64> {
    let mut mean = vec![0.0; arena.get(start).composition.len()];
    for index in start..end {
        for (sum, value) in mean.iter_mut().zip(&arena.get(index).composition) {
            *sum += value;
        }
    }
    for value in &mut mean {
        *value /= (end - start) as f64;
    }
    let total = mean.iter().sum::<f64>();
    for value in &mut mean {
        *value /= total;
    }
    mean
}

fn jensen_shannon(left: &[f64], right: &[f64]) -> f64 {
    0.5 * left
        .iter()
        .zip(right)
        .map(|(left, right)| {
            let middle = 0.5 * (left + right);
            let a = if *left == 0.0 {
                0.0
            } else {
                left * (left / middle).ln()
            };
            let b = if *right == 0.0 {
                0.0
            } else {
                right * (right / middle).ln()
            };
            a + b
        })
        .sum::<f64>()
}

const fn valid_nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
const fn valid_positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum TrajectoryObserverError {
    #[error("invalid trajectory-observation policy: {0}")]
    InvalidPolicy(&'static str),
    #[error("iterations must increase strictly")]
    NonIncreasingIteration,
    #[error("physical time must be finite and increase strictly when present")]
    InvalidPhysicalTime,
    #[error("abundance must be nonempty, nonnegative, finite, and have positive finite mass")]
    InvalidAbundance,
    #[error("{name} dimension changed from {expected} to {actual}")]
    DimensionChanged {
        name: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error("equilibrium evidence is invalid or does not match the detector observable")]
    InvalidEvidence,
    #[error("the presence of a distinct detector observable changed during observation")]
    DetectorObservableModeChanged,
    #[error("no observations are available")]
    NoSamples,
    #[error("the termination signal does not describe the final iteration")]
    StopReasonIterationMismatch,
    #[error("the retained history does not contain the accepted detection interval")]
    DetectionHistoryUnavailable,
    #[error(transparent)]
    TerminalState(#[from] TerminalStateError),
}
