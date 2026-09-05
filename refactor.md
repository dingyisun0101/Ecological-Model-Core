# Eco Core compatibility with Workflow 0.13.5

Reviewed on 2026-09-05: local `main`, commit
`6ca9e6afb3b76d8dbbd7cd343b9f5781daacc554`. Application source was read only and
validation ran in an isolated copy.

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

## Impact and recommended update

`ecological-state-toolkit` depends on Workflow 0.13.4 to provide ecological state
schema integration. Update that declaration and the lockfile to 0.13.5 when
coordinating the downstream release. The current semver range allows 0.13.5,
but an existing lockfile can retain 0.13.4; verify the actual resolved version.

State schemas, static providers, payload ownership, and ordinary observation
interfaces are unchanged. Eco Core has no current raw dependency-context parser
or old Workflow Python import requiring migration. Do not add a Python runtime
dependency to the Rust library solely for this refactor. The environment setup
above applies when an application uses `$npy` or the Python companion.

Keep ecological schemas and conversion rules in Eco Core. New Workflow typed
selectors and project helpers belong at application orchestration boundaries;
there is no reason to thread those mechanisms through schema providers.

## Validation evidence and limits

An isolated copy with the qualified Workflow 0.13.5 source override passed
**all 32 Rust tests**, without changing Eco Core source. Repeat
`cargo test --workspace --all-targets` with the published version selected in the
lockfile. GLV and Simulator's own lockfiles may select a different Eco Core patch;
check the coordinated application graph separately. This task adds only this
migration/compatibility guide and does not publish a new Eco Core version.
