# Eco Core migration to Workflow 0.13.5 and PiP 4.1.0-alpha

Initial audit on 2026-09-05: local `main`, commit
`6ca9e6afb3b76d8dbbd7cd343b9f5781daacc554`, reviewed without source edits and
validated in an isolated copy. The subsequent 0.13.2 migration and validation
are recorded below.

## Coordinated release requirements

[Workflow v0.13.5](https://github.com/dingyisun0101/Scientific-Workflow/releases/tag/v0.13.5)
publishes Rust `scientific-workflow` **0.13.5** and Python companion
`scientific-workflow` **0.4.3** (import `scientific_workflow`). Macros stay 0.2.1.
**The requested patch release contains breaking changes; no old Python aliases exist.**

**LINUX ONLY. Workflow's Python utilities require Python 3.14+. Upgrade older
Python versions. ACTIVATE THE ENVIRONMENT BEFORE EVERY LAUNCH, INCLUDING EACH
NEW SHELL. Cargo does not install or activate Python.**

```sh
python3.14 -m venv .venv
source .venv/bin/activate
python -m pip install \
  'scientific-workflow[npy] @ git+https://github.com/dingyisun0101/Scientific-Workflow.git@v0.13.5#subdirectory=python'
```

Omit `[npy]` for core recording readers without conversion/NPY views. Python is
published through the Git tag and GitHub wheel/source assets, not claimed on PyPI.
Rust 0.13.5 is published on crates.io.

**REQUIRED LAYOUT: `<study>/wf_configs/study.json` and `parameters.json`.
DO NOT RENAME OR RELOCATE THEM.** Accessors assume the standard layout. Programs
should read resolved parameter snapshots through the supported accessors, so
sweeps and overrides are respected.

Periodic recordings remain format 7; explicit `initial_and_final` recordings use
format 8. Both new readers accept 7 and 8. **NPY remains v2** and project schema
remains 1. Read the [upstream migration guide](https://github.com/dingyisun0101/Scientific-Workflow/blob/v0.13.5/docs/migration-0.13.5.md)
and [API references](https://github.com/dingyisun0101/Scientific-Workflow/blob/v0.13.5/README.md#subsystem-contracts).

## Implemented dependency updates

`ecological-state-toolkit` 0.13.2 updates its Workflow dependency declaration and
lockfile from 0.13.4 to 0.13.5 to provide ecological state schema integration.
Downstream application lockfiles must also be refreshed; verify the actual
resolved version rather than relying only on the manifest's semver range.

The release also pins `physics_in_parallel = "=4.1.0-alpha"`, following the
upstream PiP release and yanking of earlier versions. Existing lockfiles can
still use the yanked 4.0.0-alpha.2, but fresh resolution cannot. Downstream
applications exchanging PiP types with Eco Core must migrate to the same PiP
version. PiP retains the schema-v2 tensor wire format; its stricter validation
and numerical changes are described in the
[PiP migration notes](https://github.com/dingyisun0101/Physics-in-Parallel/blob/v4.1.0-alpha/docs/RELEASES.md).

State schemas, static providers, payload ownership, and ordinary observation
interfaces are unchanged. Eco Core has no current raw dependency-context parser
or old Workflow Python import requiring migration. Do not add a Python runtime
dependency to the Rust library solely for this refactor. The environment setup
above applies when an application uses `$npy` or the Python companion.

Keep ecological schemas and conversion rules in Eco Core. New Workflow typed
selectors and project helpers belong at application orchestration boundaries;
there is no reason to thread those mechanisms through schema providers.

## Validation evidence and limits

The initial isolated audit with the qualified Workflow 0.13.5 source override
passed **all 32 Rust tests**, without changing Eco Core source.

The 0.13.2 migration on 2026-09-05 also passed **all 32 Rust tests** with
`cargo test --locked --workspace --all-targets` against the published crates.io
Workflow 0.13.5 and PiP 4.1.0-alpha packages, without source overrides. Dependency
tree checks confirm one version of each. `cargo test --locked --doc` passed
four documentation tests, with one illustrative example ignored. Formatting
validation passed with `cargo fmt --all -- --check`.

`cargo publish --dry-run --locked --allow-dirty` verified the release package.
`cargo publish --locked --allow-dirty` then published `ecological-state-toolkit`
0.13.2 to crates.io and confirmed registry availability on 2026-09-05.

GLV and Simulator's own lockfiles may select a different Eco Core patch; check
the coordinated application graph separately when migrating those applications.
