# Ecological State Toolkit Design

## Status and scope

This document defines the clean-slate public `ecological-state-toolkit` crate.
There is no compatibility requirement for the former
`ecological-initial-state` crate, historical terminal-state documents, or old
consumer recordings.

The repository owns neutral ecological inputs, products, trajectory analysis,
and generic conversion of completed ecological recordings. It does not own a
model, a Workflow runtime, a study's recording selection policy, or
application orchestration.

The package and Rust library names are:

```toml
name = "ecological-state-toolkit"
```

```rust
ecological_state_toolkit
```

## Public modules

```text
ecological_state_toolkit
├── artifact
├── initial_state
├── inputs
├── interaction
├── state_schema
├── trajectory
└── terminal_state
```

The root should expose modules rather than flattening all domain types into a
large prelude.

### `artifact`

This module exposes the shared immutable descriptor, `Created`/`Reused`
publication result, and typed publication/load errors. Generic byte
publication and loading remain crate-private: `initial_state` and
`interaction` are the only public semantic write/read boundaries.

Artifacts are rooted at an explicit caller-owned `&Path`, stored beneath its
`inputs/` child, published atomically under a SHA-256-derived filename, and
loaded only after normalized-relative-path, canonical containment, and digest
verification. This facility owns prepared ecological inputs, never Workflow
model recordings.

### `initial_state`

This module owns categorical ecological initialization:

- categorical taxon distributions;
- random and centered-seed generation;
- validated taxon counts and lattice state;
- generated-state RNG provenance;
- a single source contract for generated or verified inputs;
- stable initial-state JSON and content-addressed ecological input artifacts;
- verified loading and persistence.

Each model owns conversion from a categorical state into its native evolving
state. No model-specific conversion belongs here.

### `interaction`

This module owns model-neutral ecological interactions:

- the interaction matrix backed by PiP dense matrix storage;
- matrix dimension and finite-value validation;
- generated-matrix recipes and RNG provenance;
- PiP-backed transpose, subtraction, scaling, and maximum-magnitude reduction;
- model-neutral derived transformations with complete provenance;
- stable interaction JSON;
- content-addressed ecological input artifacts;
- verified loading and persistence.

The primitive independent recipes sample every matrix entry, including the
diagonal, without hidden normalization or sparsity policy:

```text
RandomUniform:  A[i,j] ~ Uniform(minimum, maximum)
RandomGaussian: A[i,j] ~ Normal(mean, standard_deviation)
```

Model preparation is expressed as explicit, non-mutating transformations.
`antisymmetrize()` returns exactly `A - A^T`, `scale(scalar)` returns
`scalar * A`, and `normalize(threshold)` scales down only when the maximum
absolute entry exceeds the threshold. These operations compose PiP matrix
APIs; this crate does not reimplement transposition, elementwise subtraction,
scalar multiplication, or matrix reduction.

Recipes, rather than `Config`-named domain types, describe scientific
construction. The public vocabulary is `InitialStateRecipe`,
`InteractionMatrixRecipe`, policy, descriptor, and provenance. `Config`
remains reserved for upstream configuration types such as `RngConfig` and
`SquareLatticeConfig`.

Specialized correlated and sign-structured ecological ensembles remain
available and keep their explicit normalization, diagonal, reciprocal-pair,
connectance, and randomness policies. Applications treat recipes opaquely:
they request a source realization and compose transformations without matching
a specific recipe family.

Each consumer owns its model-specific interpretation, such as direct
continuous-time coefficients or conversion into replacement probabilities.
The shared matrix must not impose either contract. Probability clamping is a
model transformation, not an ecological matrix-generation recipe.

### `trajectory`

This module owns bounded, model-independent observation of an ecological
trajectory. It receives borrowed values from a model after completed steps and
may signal a numerically accepted long-term behavior. It never owns or reads a
model state, advances time, mutates a runtime context, or completes a recording.

The central type is `TrajectoryObserver`. `StateObserver` is intentionally not
used because that term commonly means hidden-state reconstruction in control
theory.

### `terminal_state`

This module owns the terminal scientific product constructed by
`TrajectoryObserver::finish`:

- normalized terminal composition;
- final iteration and optional physical time;
- termination classification and reason;
- represented sample interval and count;
- equilibrium, periodic-orbit, and absorption diagnostics when applicable;
- strict serialization and validation.

This is a new format. No legacy parser or format alias is required.

## Dependency direction

```text
Physics in Parallel ─┐
Serde / SHA-256 ─────┼──> Rust crate ───────────────> model crates
Scientific Workflow ┘          │
                               │ shared ecological contracts
Workflow Reader + NumPy ─────> Python companion ───> study adapters / analysis
```

The Rust crate must not depend on a downstream model, application
orchestrator, ndarray, or another runtime implementation. Its only Workflow
dependency is the neutral `StateSchemaProvider` type used to hand the embedded
schema to a downstream execution unit.

The Rust crate must not decode, resolve, or inspect Workflow task
configuration. Its public scientific APIs accept resolved in-memory values.
Direct data and artifact I/O accepts ordinary paths or descriptors;
configuration keys are resolved by an orchestrator before the core is called.

Complete scientific applications deliberately sit on the other side of this
boundary and own canonical workload directories, including model inputs,
sweeps, path resolution, and writer settings. They translate those documents
into the plain core recipes, observations, and artifacts defined here.

The artifact, initial-state, and interaction modules use only ordinary paths,
toolkit descriptors, and fully resolved PiP `RngConfig` values. No
Workflow-owned type crosses their public API. The trajectory and
terminal-state modules depend only on plain scientific values and time
coordinates.

The `inputs` module composes, but does not reinterpret, two existing artifact
contracts. `EcologicalInputs` contains one final model-ready interaction
reference and one canonical initial-state reference. It validates reference
envelopes and taxon dimensionality, then resolves both through their owning
modules. An application may reuse the same initial-state reference across
models while supplying different explicitly derived interaction references.
This composition does not make Ecological State Toolkit an orchestrator and
introduces no Workflow configuration, scheduling, or recording dependency.

Generated values retain the resolved PiP method and seed rather than wrapping
them in a second generic RNG record. Generator identity, version, and recipe
remain in the toolkit's domain provenance. This is sufficient for
deterministic replay without making the toolkit a metadata framework.

## Python recording conversion

The installable `python/ecological_state_toolkit` companion is the generic
analysis conversion boundary. It may depend on NumPy and the public Scientific
Workflow reader because it operates after a recording is finalized; those
dependencies do not cross into the Rust model crate.

Its public request vocabulary is `RecordingSpec -> StreamSpec -> FieldSpec`.
Fields select stable ecological encodings, never a downstream execution-unit
name. It owns chunk verification through the official reader, bounded
multiprocessing, shape/dtype invariants, NPY memmap allocation, atomic member
publication, checksummed resume, descriptors, and the generic CLI request and
manifest formats.

A downstream study owns which recordings and streams to select, how a model's
field maps to a generic ecological encoding, and all scientific joining or
pairing metadata. It may retain the toolkit's result directly or publish a
study-specific manifest. The toolkit must contain no Dispatcher, GLV,
Simulator, task identity, parameter ordinal, or study-format constant.

The Python CLI is a thin transport over the library. Conversion itself is not
a Workflow execution unit: finalized-recording transformation has no evolving
state or `step()` lifecycle. Scheduling belongs to an optional downstream
program/task adapter, which invokes the unchanged library API.

## Trajectory observation modes

The common policy has three explicit modes:

```rust,ignore
pub enum TrajectoryObservationPolicy {
    Disabled,
    TerminalOnly(TerminalPolicy),
    Detect(DetectionPolicy),
}
```

`Disabled` means complete absence of the facility:

- no observer construction;
- no history or scratch allocation;
- no normalization;
- no residual calculation;
- no terminal-state product or observer metadata.

`TerminalOnly` retains only the bounded history required to construct a
terminal composition.

`Detect` retains terminal history and enables configured equilibrium,
periodic-orbit, and absorption detectors.

The core crate defines no model-wide default. A deterministic model may select
`Detect`, while a stochastic or discrete model may prefer `TerminalOnly` or
`Disabled`. Consumers derive any default observation interval from their own
recording policy and pass it to the core as a plain value; the core does not
inspect Workflow stream configuration. An explicit observer policy may
override the interval. Final observation is always included by `finish`, even
when it is off cadence.

## Borrowed observation contract

The engine calls `observe` after each successfully completed scientific step.
The input is borrowed for that call:

```rust,ignore
pub struct TrajectoryObservation<'a> {
    pub iteration: u64,
    pub physical_time: Option<f64>,
    pub abundance: AbundanceView<'a>,
    pub detector_observable: Option<&'a [f64]>,
    pub equilibrium_evidence: EquilibriumEvidence<'a>,
}

pub enum AbundanceView<'a> {
    Continuous(&'a [f64]),
    Counts(&'a [usize]),
}
```

Continuous and count inputs are validated and normalized into the same
internal composition. `Continuous` covers relative frequencies, densities,
and continuous absolute populations; `Counts` preserves an exact integer input
boundary. Original mass is retained as a scalar for mass-stability checks.

`detector_observable` defaults to normalized abundance. GLV may instead supply
a borrowed flattened spatial field, preventing a spatially evolving state with
constant aggregate abundance from being classified as stationary. Terminal
composition always comes from `abundance`, never from the optional detector
field.

An observation must have increasing iteration, stable abundance dimension,
valid time, nonempty nonnegative finite values, and positive finite mass.
Detector dimension must remain stable while its configured history is active.

## Equilibrium evidence

A quiet sampled trajectory is insufficient evidence of an equilibrium. The
observer accepts an equilibrium only when window evidence and
model-authoritative evidence both pass.

```rust,ignore
pub enum EquilibriumEvidence<'a> {
    Unavailable,
    Residual { values: &'a [f64] },
    MaximumScaledResidual { value: f64 },
    AbsorbingState,
}
```

`Residual` is preferred for ODE/PDE models. It is the borrowed authoritative
vector field or discretized PDE residual for the exact submitted detector
observable. The observer owns the common component-wise scaling:

```text
max_i |r_i| / (absolute_tolerance + relative_tolerance * |x_i|)
```

The residual and observable must have equal length and all residual components
must be finite.

`MaximumScaledResidual` permits a kernel to compute the same maximum without
exposing or allocating a residual vector. Its value must be for the exact
submitted state, use the configured tolerances, cover every authoritative
component, and be finite and nonnegative.

`AbsorbingState` is a model assertion that no future transition can leave the
submitted state. The observer validates its abundance input but the model owns
the no-transition proof. It is not treated as interchangeable with a small
continuous vector-field residual.

`Unavailable` prevents equilibrium acceptance whenever the equilibrium
detector requires evidence at that sample. Evidence is ignored on steps that
the observer does not retain. A consumer should query cadence/evidence demand
before performing an expensive residual evaluation.

## Long-term behavior terminology

Public terminology follows dynamical-systems and theoretical-ecology usage:

- `Equilibrium`: a steady state with vanishing equations of motion;
- `PeriodicOrbit`: a numerically recurrent closed trajectory;
- `AbsorbingState`: a discrete or stochastic state that cannot be left;
- `TrailingAverage`: a bounded terminal estimate without an accepted
  asymptotic classification.

The crate must not call a detected periodic orbit a `LimitCycle`, because the
detector does not establish isolation or attraction. It must not call it a
`NeutralCycle`, because recurrence does not establish neutral stability.

Signals are model-independent and contain auditable diagnostics:

```rust,ignore
pub enum TerminationSignal {
    Equilibrium(EquilibriumDiagnostics),
    PeriodicOrbit(PeriodicOrbitDiagnostics),
    AbsorbingState(AbsorptionDiagnostics),
}
```

These are numerical acceptances under declared tolerances, not proofs of
linear, nonlinear, ecological, or evolutionary stability.

## Detection semantics

### Equilibrium

Equilibrium acceptance retains the current GLV scientific requirements:

- invariant support throughout each confirmation window;
- bounded composition variation;
- optional bounded relative mass variation;
- authoritative scaled residual at or below one;
- staged, increasing confirmation windows;
- exact absorbing shortcut only with `AbsorbingState` evidence.

A support change resets staged equilibrium confirmation.

### Periodic orbit

Periodic-orbit acceptance retains the current recurrence requirements:

- period within configured sample bounds;
- at least the configured number of repeated cycles;
- bounded recurrence distance between corresponding cycle samples;
- nontrivial minimum cycle amplitude.

The result records period in samples, represented iteration span, repeated
cycles, maximum recurrence distance, and amplitude. It does not claim orbit
stability. `finish` uses complete detected cycles to construct the orbit mean.

### Absorption

Absorption requires explicit `AbsorbingState` evidence. Support size alone is
not a generic proof because different models may permit invasion, mutation,
forcing, or other transitions from a single-supported state.

## Finalization

Every successful observed run calls `finish`, regardless of why evolution
stopped:

```rust,ignore
observer.finish(final_observation, StopReason::Detected(signal))
observer.finish(final_observation, StopReason::MaximumIterations)
observer.finish(final_observation, StopReason::Requested)
observer.finish(final_observation, StopReason::ModelSpecific(reason))
```

Finalization consumes the observer and returns one validated `TerminalState`.
It forces the final observation into history if it is not already present.

- equilibrium: exact final normalized composition;
- absorbing state: exact final normalized composition;
- periodic orbit: mean over complete detected cycles;
- maximum/requested/model-specific stop: bounded trailing mean.

In `Disabled` mode there is no observer and therefore no terminal product.
Consumers may still record their own stop reason and completed iteration.

## Allocation and cloning contract

Calling `observe` on an off-cadence step performs validation only as needed for
control flow and owns no model data. `Disabled` performs no observer work at
all.

At retained samples, copying is unavoidable because the model mutates its
buffers on the next step. The implementation minimizes it as follows:

- use a fixed-capacity circular arena of reusable sample slots;
- allocate slot vectors after the first retained observation establishes
  dimensions;
- normalize directly into a slot without an intermediate vector;
- reuse slot buffers on overwrite;
- retain one normalized abundance composition per sample;
- alias detector logic to that composition when no separate observable exists;
- allocate a second retained vector only for a distinct detector observable;
- store mass, time, evidence score, and sequence position as scalars;
- store support in reusable compact masks or derive it without duplicating the
  full sample;
- let all detectors and the terminal estimator address the same arena;
- reuse one final accumulator and move it into `TerminalState`.

The observer never clones `SystemState`, ndarray arrays, PiP lattices, or a
complete simulation. The expected ownership path is one copy per retained
sample and one final vector moved into the terminal product.

## Consumer/runtime boundary

The observer returns a signal only. The consuming model or host decides whether
to stop, updates its runtime progress, completes recording, and publishes
products. Core errors contain scientific context but no task or phase status.

## Validation strategy

The crate requires independent tests for:

- frequency/count normalization equivalence;
- cadence and forced-final sampling;
- zero observer allocations in `Disabled` construction paths;
- stable retained allocation capacity after warm-up;
- no duplicate detector vector for abundance-based detection;
- equilibrium rejection without residual evidence;
- raw and pre-scaled residual equivalence;
- support-change reset;
- staged equilibrium acceptance;
- periodic-orbit period, recurrence, and amplitude diagnostics;
- rejection of constant trajectories as periodic orbits;
- absorption only with explicit evidence;
- terminal exact-state, orbit-mean, and trailing-mean construction;
- spatial detector observables independent from terminal abundance;
- strict serialization and invalid-product rejection;
- initial-state and interaction artifact round trips;
- atomic reuse, corruption rejection, and artifact-root containment;
- deterministic resolved RNG provenance.

Downstream model crates retain end-to-end equivalence, recording, and runtime
tests.
