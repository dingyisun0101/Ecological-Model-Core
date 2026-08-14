//! Validated terminal ecological composition and long-term-behavior diagnostics.

use serde::{Deserialize, Serialize};
use thiserror::Error;

pub const TERMINAL_STATE_FORMAT: &str = "ecological.terminal-state.v1";
pub const TERMINAL_STATE_METADATA_KEY: &str = "terminal_state";

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct EquilibriumDiagnostics {
    pub iteration: u64,
    pub completed_windows: usize,
    pub final_window_samples: usize,
    pub maximum_observable_distance: f64,
    pub relative_mass_range: f64,
    pub maximum_scaled_residual: f64,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct PeriodicOrbitDiagnostics {
    pub iteration: u64,
    pub period_samples: usize,
    pub repeated_cycles: usize,
    pub first_cycle_iteration: u64,
    pub last_cycle_iteration: u64,
    pub maximum_recurrence_distance: f64,
    pub orbit_amplitude: f64,
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct AbsorptionDiagnostics {
    pub iteration: u64,
    pub supported_taxa: usize,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "kind", content = "diagnostics", rename_all = "snake_case")]
pub enum TerminationSignal {
    Equilibrium(EquilibriumDiagnostics),
    PeriodicOrbit(PeriodicOrbitDiagnostics),
    AbsorbingState(AbsorptionDiagnostics),
}

impl TerminationSignal {
    pub const fn iteration(&self) -> u64 {
        match self {
            Self::Equilibrium(value) => value.iteration,
            Self::PeriodicOrbit(value) => value.iteration,
            Self::AbsorbingState(value) => value.iteration,
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(tag = "kind", content = "value", rename_all = "snake_case")]
pub enum StopReason {
    Detected(TerminationSignal),
    MaximumIterations,
    Requested,
    ModelSpecific(String),
}

impl StopReason {
    pub const fn signal(&self) -> Option<&TerminationSignal> {
        match self {
            Self::Detected(signal) => Some(signal),
            _ => None,
        }
    }
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TerminalClassification {
    Equilibrium,
    PeriodicOrbit,
    AbsorbingState,
    TrailingAverage,
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct TerminalState {
    format: String,
    classification: TerminalClassification,
    stop_reason: StopReason,
    iteration: u64,
    physical_time: Option<f64>,
    composition: Vec<f64>,
    sample_count: usize,
    first_sample_iteration: u64,
    last_sample_iteration: u64,
}

impl TerminalState {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        classification: TerminalClassification,
        stop_reason: StopReason,
        iteration: u64,
        physical_time: Option<f64>,
        composition: Vec<f64>,
        sample_count: usize,
        first_sample_iteration: u64,
        last_sample_iteration: u64,
    ) -> Result<Self, TerminalStateError> {
        let state = Self {
            format: TERMINAL_STATE_FORMAT.to_owned(),
            classification,
            stop_reason,
            iteration,
            physical_time,
            composition,
            sample_count,
            first_sample_iteration,
            last_sample_iteration,
        };
        state.validate()?;
        Ok(state)
    }

    pub fn format(&self) -> &str {
        &self.format
    }
    pub const fn classification(&self) -> TerminalClassification {
        self.classification
    }
    pub const fn stop_reason(&self) -> &StopReason {
        &self.stop_reason
    }
    pub const fn iteration(&self) -> u64 {
        self.iteration
    }
    pub const fn physical_time(&self) -> Option<f64> {
        self.physical_time
    }
    pub fn composition(&self) -> &[f64] {
        &self.composition
    }
    pub const fn sample_count(&self) -> usize {
        self.sample_count
    }
    pub const fn first_sample_iteration(&self) -> u64 {
        self.first_sample_iteration
    }
    pub const fn last_sample_iteration(&self) -> u64 {
        self.last_sample_iteration
    }
    pub fn to_json_bytes(&self) -> Result<Vec<u8>, TerminalStateError> {
        Ok(serde_json::to_vec(self)?)
    }
    pub fn from_json_bytes(bytes: &[u8]) -> Result<Self, TerminalStateError> {
        let state: Self = serde_json::from_slice(bytes)?;
        state.validate()?;
        Ok(state)
    }

    pub fn validate(&self) -> Result<(), TerminalStateError> {
        if self.format != TERMINAL_STATE_FORMAT
            || self.sample_count == 0
            || self.first_sample_iteration > self.last_sample_iteration
            || self.last_sample_iteration > self.iteration
            || self.physical_time.is_some_and(|value| !value.is_finite())
        {
            return Err(TerminalStateError::InvalidProduct);
        }
        validate_composition(&self.composition)?;
        match (&self.classification, &self.stop_reason) {
            (
                TerminalClassification::Equilibrium,
                StopReason::Detected(TerminationSignal::Equilibrium(value)),
            ) if value.iteration == self.iteration && self.sample_count == 1 => {}
            (
                TerminalClassification::PeriodicOrbit,
                StopReason::Detected(TerminationSignal::PeriodicOrbit(value)),
            ) if value.iteration == self.iteration
                && self.sample_count >= value.period_samples
                && self.first_sample_iteration == value.first_cycle_iteration
                && self.last_sample_iteration == value.last_cycle_iteration => {}
            (
                TerminalClassification::AbsorbingState,
                StopReason::Detected(TerminationSignal::AbsorbingState(value)),
            ) if value.iteration == self.iteration && self.sample_count == 1 => {}
            (
                TerminalClassification::TrailingAverage,
                StopReason::MaximumIterations | StopReason::Requested,
            ) => {}
            (TerminalClassification::TrailingAverage, StopReason::ModelSpecific(value))
                if !value.trim().is_empty() => {}
            _ => return Err(TerminalStateError::ClassificationMismatch),
        }
        Ok(())
    }
}

fn validate_composition(values: &[f64]) -> Result<(), TerminalStateError> {
    if values.is_empty()
        || values
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(TerminalStateError::InvalidComposition);
    }
    let total = values.iter().sum::<f64>();
    if !total.is_finite() || (total - 1.0).abs() > 1.0e-10 {
        return Err(TerminalStateError::InvalidComposition);
    }
    Ok(())
}

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum TerminalStateError {
    #[error("terminal-state document is structurally invalid")]
    InvalidProduct,
    #[error("terminal composition must be nonempty, finite, nonnegative, and normalized")]
    InvalidComposition,
    #[error("terminal classification and stop reason are inconsistent")]
    ClassificationMismatch,
    #[error(transparent)]
    Json(#[from] serde_json::Error),
}
