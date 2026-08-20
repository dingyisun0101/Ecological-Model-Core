# Ecological Model Core

Shared scientific primitives for ecological models that use different numerical
methods. The crate deliberately contains no simulation engine or Workflow task:
GLV and lattice Simulator remain independent runtime-integrated model crates.

The public modules are:

- `initial_state`: reproducible categorical lattice initial states and verified
  Workflow artifacts;
- `interaction`: validated ecological interaction matrices, reproducible random
  sources, and composable PiP-backed transformations;
- `trajectory`: a bounded, allocation-conscious trajectory observer with
  disabled, terminal-only, and detection modes;
- `terminal_state`: validated terminal composition and auditable equilibrium,
  periodic-orbit, and absorbing-state results.

Construction parameters are called *recipes*, keeping scientific realization
instructions distinct from Workflow configuration. Engines submit borrowed
abundance slices, so continuous and lattice models retain their native state
representations and only copy the bounded samples the selected observer mode
needs.

The crate never reads or resolves Workflow task configuration. Callers provide
ordinary resolved values (`Vec`, slices, `PathBuf`, recipes, and descriptors);
an orchestrator or example owns configuration decoding and path-key resolution.

Interaction matrices are sampled before model-specific transformations are
applied. Every entry of `RandomUniform` and `RandomGaussian`, including the
diagonal, is sampled independently. Transformations return new matrices, leave
their sources unchanged, and retain a complete derived-provenance chain:

```rust
use ecological_model_core::interaction::{InteractionMatrix, InteractionMatrixRecipe};
use physics_in_parallel::rng::RngConfig;

let recipe = InteractionMatrixRecipe::RandomGaussian {
    mean: 0.0,
    standard_deviation: 1.0,
    rng: RngConfig::new(Some(42), None),
};
let lattice = InteractionMatrix::generate(8, &recipe)?
    .scale(0.5)?
    .abs()?;
lattice.ensure_max_abs_at_most(1.0)?;
let glv = lattice.antisymmetrize()?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

`abs` applies elementwise absolute value. `clamp_min` and `clamp_max` provide
independent finite lower and upper bounds. `ensure_max_abs_at_most` validates
a bound without changing the matrix, which is useful when exceeding the
threshold should be an error.

```rust
use ecological_model_core::trajectory::{
    TerminalPolicy, TrajectoryObservationPolicy, TrajectoryObserver,
};

let observer = TrajectoryObserver::from_policy(
    TrajectoryObservationPolicy::TerminalOnly(TerminalPolicy {
        sample_interval_iterations: 100,
        trailing_window_samples: 20,
    }),
)?;
assert!(observer.is_some());
# Ok::<(), Box<dyn std::error::Error>>(())
```

The crate is a clean-start API. It does not read artifacts produced by the
former `ecological-initial-state` crate.
