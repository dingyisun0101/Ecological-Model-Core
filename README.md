# Ecological Model Core

> **Alpha and breaking API notice:** version 0.10.0 targets PiP 3.7 and Workflow
> 0.9. Earlier `eco_core` releases use obsolete dependency contracts and are
> unsupported.

Shared scientific primitives for ecological models that use different numerical
methods. The crate deliberately contains no simulation engine or Workflow task:
each downstream model remains an independent runtime-integrated crate.

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

`InitialStateRecipe::BalancedUniform` constructs the unique minimum-spread
uniform count allocation: every taxon receives either `floor(N/K)` or
`ceil(N/K)` sites. PiP then applies one reproducible unbiased permutation to
the complete row-major lattice tensor. The resulting counts differ by at most
one, and `InitialState::frequencies` exposes their exact aggregate composition
for a continuous model using the same initial condition.

The crate never reads or resolves Workflow task configuration. Callers provide
ordinary resolved values (`Vec`, slices, `PathBuf`, recipes, and descriptors);
an orchestrator or example owns configuration decoding and path-key resolution.

Interaction matrices are sampled before model-specific transformations are
applied. Every entry of `RandomUniform` and `RandomGaussian`, including the
diagonal, is sampled independently. Transformations return new matrices, leave
their sources unchanged, and retain a complete derived-provenance chain:

```rust
use ecological_model_core::interaction::InteractionMatrixRecipe;
use physics_in_parallel::prelude::basic::RngConfig;

let recipe = InteractionMatrixRecipe::RandomGaussian {
    mean: 0.0,
    standard_deviation: 1.0,
    rng: RngConfig::new(Some(42), None),
};
let lattice = recipe.generate(8)?
    .scale(0.5)?
    .abs()?
    .normalize(1.0)?;
lattice.ensure_max_abs_at_most(1.0)?;
let glv = lattice.antisymmetrize()?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

`abs` applies elementwise absolute value. `clamp_min` and `clamp_max` provide
independent finite lower and upper bounds. `ensure_max_abs_at_most` validates
a bound without changing the matrix, which is useful when exceeding the
threshold should be an error.

Matrix constructors infer the species count from the square PiP matrix; callers
do not repeat it. Only random recipes require a species count because it defines
the matrix that has not yet been created. `InteractionMatrix` delegates scalar,
transpose, subtraction, absolute-value, reduction, and vector-application work
to PiP. `mul_vectors_into` applies the same matrix to a contiguous batch while
reusing caller-owned output storage.

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
