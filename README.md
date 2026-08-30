# Ecological State Toolkit

> **Breaking 0.12 update:** `ecological-model-core` is superseded by the
> `ecological-state-toolkit` crate and `ecological_state_toolkit` Rust/Python
> import. There are no compatibility aliases. The ecological schema provider
> identity is now `ecological-state-toolkit.ecological-state.v1`.

An end-to-end, model-neutral package for ecological state construction,
validation, schema provision, trajectory products, and conversion of completed
Workflow recordings into analysis-efficient NumPy arrays. It deliberately
contains no simulation engine or Workflow task: each downstream model and study
remains independently runtime-integrated.

The public modules are:

- `artifact`: shared immutable descriptor, publication disposition, and typed
  integrity failures; byte publication itself remains behind the ecological
  modules that know the document semantics;
- `initial_state`: reproducible categorical lattice initial states and verified
  ecological input artifacts;
- `inputs`: the canonical validated `(InteractionMatrix, InitialState)`
  envelope shared by independent engines;
- `interaction`: validated ecological interaction matrices, reproducible random
  sources, and composable PiP-backed transformations;
- `state_schema`: the sole embedded canonical ecological Workflow state schema
  and its typed provider descriptor;
- `trajectory`: a bounded, allocation-conscious trajectory observer with
  disabled, terminal-only, and detection modes;
- `terminal_state`: validated terminal composition and auditable equilibrium,
  periodic-orbit, and absorbing-state results.

The companion Python module `ecological_state_toolkit` exposes generic,
verified recording conversion; its API is documented below.

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

The Rust crate depends on Scientific Workflow only for the neutral
`StateSchemaProvider` descriptor. It never reads or resolves a Workflow
project, implements an execution unit, or writes a Workflow recording. It has
no dependency on Dispatcher, GLV, Simulator, or another private downstream
crate. Callers provide ordinary resolved values (`Vec`, slices, `PathBuf`,
recipes, and descriptors); an orchestrator or model owns configuration
decoding and path-key resolution.

For comparative applications, `inputs::EcologicalInputs` is the canonical
model-neutral envelope. It pairs one final model-ready
`InteractionArtifactReference` with one `InitialStateArtifactReference`,
validates their shared taxon dimension without IO, and resolves both with full
digest verification. Multiple models can receive the same initial-state
reference: a categorical model consumes its space and counts, while a
continuous model derives exact frequencies from that same state. The envelope
does not select recipes, apply model-specific transformations, or own runtime
configuration.

After verification, `ResolvedEcologicalInputs` is exactly the typed wrapper
around `(InteractionMatrix, InitialState)`. Models may borrow the two values or
call `into_parts()` to take that tuple directly; Ecological State Toolkit does
not construct a GLV- or Simulator-specific state around it.

## Installation

```toml
[dependencies]
ecological-state-toolkit = "0.12.2"
scientific-workflow = "0.11.5"
physics_in_parallel = "3.7.0"
```

The companion Python distribution is installed from the same tagged source:

```text
python -m pip install \
  "ecological-state-toolkit @ git+https://github.com/dingyisun0101/Ecological-Model-Core.git@v0.12.2#subdirectory=python"
```

Use this crate when multiple ecological models need the same validated
initial-state, interaction, trajectory, or terminal-product semantics. A model
that needs only one small local calculation may be clearer without the extra
dependency.

When used with Scientific Workflow 0.11.5, put recipes and other resolved
scientific values in the model's custom `Constants` type. The registered model
still directly owns its Workflow `SystemState`; Ecological State Toolkit owns
only the standard layout supplied to that state, not the model, observation
plan, recording destination, or runtime.

## Recording conversion

The Python distribution is the analysis-facing half of the toolkit. It uses
the official Workflow reader to verify every JSONL chunk and writes each
selected stream directly into atomic, C-contiguous NPY memmaps. Conversion is
bounded by a caller-selected process count and resumes at recording
granularity only when both the source metadata checksum and complete
conversion request still match.

Applications use the library API as the stable boundary:

```python
from pathlib import Path

from ecological_state_toolkit import (
    ArrayEncoding,
    FieldSpec,
    RecordingSpec,
    StreamSpec,
    convert_recordings,
)

specification = RecordingSpec(
    recording=Path("recording/member-000000"),
    identity="reference",
    streams=(
        StreamSpec(
            "signal",
            (
                FieldSpec("abundance", ArrayEncoding.TENSOR_F64, "values"),
                FieldSpec("total", ArrayEncoding.FLOAT_SCALAR, "total"),
            ),
        ),
    ),
    metadata={"role": "reference"},
)
converted = convert_recordings(
    (specification,), Path("processed"), workers=4
)
signal = converted[0].arrays["signal_values"]
```

`RecordingSpec` identifies one completed recording, its generic stream
contracts, and optional JSON metadata retained for its caller. `StreamSpec`
names one Workflow stream. `FieldSpec` maps a field to an output basename and
one `ArrayEncoding`: f64 tensor, nonnegative f64 tensor/vector, nonnegative
u32 vector, categorical lattice, finite float scalar, or nonnegative integer
scalar. Categorical
lattices additionally require `category_count`, which determines validation
and the smallest safe unsigned storage dtype. `convert_recordings` takes a
sequence of specifications, an output `Path`, a positive worker count, and an
optional parent-process progress callback. It returns typed
`ConvertedRecording` values with array descriptors, source checksums, source
user metadata, and caller metadata.

For shell use, `ecological-state-convert --request REQUEST.json --output DIR`
accepts an `ecological-state-toolkit.conversion-request.v1` document and
publishes a generic batch manifest. This CLI is deliberately a thin adapter
over the library, not a raw script API.

Conversion is not an `ExecutionUnit`: it consumes finalized recordings after
scientific execution and has no evolving state, step, completion, or ensemble
semantics. A study may call the library after its export phase and then add its
own selection, pairing, and manifest semantics. If conversion later needs to
be scheduled independently, a study-owned program/task adapter should invoke
this same library without moving study semantics into the toolkit.

## Canonical ecological state schema

Ecological State Toolkit is the sole source of
`schemas/ecological_state.json`. A receiving model does not vendor or load
another copy. It forwards the provider from its Workflow execution-unit
implementation:

```rust,ignore
use ecological_state_toolkit::state_schema::ecological_state_schema;
use scientific_workflow::prelude::*;

impl ExecutionUnit for EcologicalUnit {
    type Constants = EcologicalConstants;

    fn standard_state_schema() -> Option<scientific_workflow::state::StateSchemaProvider> {
        Some(ecological_state_schema())
    }

    // preflight, initialize, member_count, member, and step follow.
}
```

The model's `study.json` task can then omit `state`, and `paths.states` can be
absent when no project-owned schema is needed. Study validates the embedded
document once, passes the resulting `SystemStateSchema` directly to the unit,
and records `ecological-state-toolkit.ecological-state.v1` as state provenance.
An explicit project state remains a deliberate override.

## Verified ecological inputs

Initial states and interaction matrices can be published as immutable,
content-addressed inputs beneath an explicit artifact root:

```rust,no_run
use std::path::Path;

use ecological_state_toolkit::artifact::ArtifactDisposition;
use ecological_state_toolkit::initial_state::{
    InitialStateRecipe, persist_initial_state, load_verified_initial_state,
};
use physics_in_parallel::prelude::basic::{RngConfig, SquareLatticeConfig};

let initial = InitialStateRecipe::BalancedUniform {
    rng: RngConfig::new(Some(42), None),
}.create(SquareLatticeConfig::periodic(&[8, 8]), 4)?;
let root = Path::new("prepared-inputs");
let persisted = persist_initial_state(root, &initial)?;
assert_eq!(persisted.disposition(), ArtifactDisposition::Created);
let restored = load_verified_initial_state(root, persisted.descriptor())?;
assert_eq!(restored.counts(), initial.counts());
# Ok::<(), Box<dyn std::error::Error>>(())
```

Publication creates `<artifact-root>/inputs`, writes and synchronizes complete
temporary bytes, then publishes a digest-derived immutable filename. Repeating
the same publication returns `Reused`; conflicting or corrupted bytes fail.
Loading rejects malformed descriptors, path traversal, symlink escape, and
SHA-256 mismatch before semantic decoding.

These are application-prepared scientific inputs, not an alternative Workflow
recording writer. During a Workflow model run, Workflow Persistence remains the
sole owner of model recording output.

## Migrating from 0.10

- Replace `scientific_workflow::ArtifactDisposition` with
  `ecological_state_toolkit::artifact::ArtifactDisposition`.
- Pass an artifact-root `&Path` to `persist_initial_state` or
  `persist_interaction_matrix`; do not construct a Workflow `ExecutionScope`.
- Pass the same root to verified loading. Serialized references now expose
  `artifact_root()` and serialize `artifact_root`, not `execution_directory`.
- Replace initial-state `rng_record()` with `rng_config()` and interaction
  `generator_rng_record()` with `generator_rng_config()`. Both return the
  resolved PiP method and seed directly.
- Remove the old initial-state and interaction RNG namespace constants. The
  domain identity/version live in the typed ecological provenance; PiP owns RNG
  method identity and seed encoding.
- Regenerate v1 initial-state documents and portable references. Version 0.11
  deliberately provides no legacy parser or alias.

Interaction matrices are sampled before model-specific transformations are
applied. Every entry of `RandomUniform` and `RandomGaussian`, including the
diagonal, is sampled independently. Transformations return new matrices, leave
their sources unchanged, and retain a complete derived-provenance chain:

```rust
use ecological_state_toolkit::interaction::InteractionMatrixRecipe;
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
let antisymmetric = lattice.antisymmetrize()?;
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
use ecological_state_toolkit::trajectory::{
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
former `ecological-initial-state` crate or `ecological-model-core` 0.10.
Initial-state documents use `ecological.initial-state.v2`; serialized
initial-state and interaction references use their respective `*.v2` formats.
