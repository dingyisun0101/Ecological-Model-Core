//! Canonical ecological state layout for Workflow-integrated models.

use scientific_workflow::state::StateSchemaProvider;

/// Stable provenance identity for the canonical ecological state layout.
pub const ECOLOGICAL_STATE_SCHEMA_ID: &str = "ecological-model-core.ecological-state.v1";

/// Returns Eco Core's canonical ecological state-schema provider.
///
/// The JSON is embedded in this crate, making Eco Core its sole owner. This
/// function performs no filesystem IO or parsing. A downstream execution unit
/// returns the descriptor from `ExecutionUnit::standard_state_schema`; Workflow
/// validates and caches it during Study preflight.
pub const fn ecological_state_schema() -> StateSchemaProvider {
    StateSchemaProvider::new(
        ECOLOGICAL_STATE_SCHEMA_ID,
        include_bytes!("../schemas/ecological_state.json"),
    )
}
