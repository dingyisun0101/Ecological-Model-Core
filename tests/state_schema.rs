use ecological_model_core::state_schema::{ECOLOGICAL_STATE_SCHEMA_ID, ecological_state_schema};

#[test]
fn canonical_ecological_schema_is_one_embedded_provider() {
    let provider = ecological_state_schema();
    assert_eq!(provider.id(), ECOLOGICAL_STATE_SCHEMA_ID);

    let document: serde_json::Value = serde_json::from_slice(provider.document()).unwrap();
    let names = document["fields"]
        .as_array()
        .unwrap()
        .iter()
        .map(|field| field["name"].as_str().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(names, ["abundance", "space", "total"]);
}
