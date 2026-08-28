//! Content-addressed ecological inputs rooted at an explicit caller-owned path.
//!
//! This module owns byte identity and integrity, not application execution or
//! Workflow persistence. Public ecological modules embed the descriptor in
//! their typed provenance and provide the semantic persistence operations.

use std::fs::{self, File, OpenOptions};
use std::io::Write as _;
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use serde::{Deserialize, Serialize};
use sha2::{Digest as _, Sha256};
use thiserror::Error;

static TEMPORARY_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// Exact identity and root-relative location of immutable ecological input.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ArtifactDescriptor {
    sha256: String,
    path: String,
}

impl ArtifactDescriptor {
    /// Returns the lowercase SHA-256 identity of the exact bytes.
    pub fn sha256(&self) -> &str {
        &self.sha256
    }

    /// Returns the normalized artifact path relative to its artifact root.
    pub fn path(&self) -> &Path {
        Path::new(&self.path)
    }
}

/// Whether publication created new bytes or reused identical existing bytes.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum ArtifactDisposition {
    /// This publication created the destination artifact.
    Created,
    /// Identical bytes already existed and were reused.
    Reused,
}

pub(crate) struct PersistedArtifact {
    descriptor: ArtifactDescriptor,
    disposition: ArtifactDisposition,
}

impl PersistedArtifact {
    pub(crate) const fn descriptor(&self) -> &ArtifactDescriptor {
        &self.descriptor
    }

    pub(crate) const fn disposition(&self) -> ArtifactDisposition {
        self.disposition
    }
}

pub(crate) struct VerifiedArtifact {
    path: PathBuf,
    bytes: Vec<u8>,
}

impl VerifiedArtifact {
    pub(crate) fn path(&self) -> &Path {
        &self.path
    }

    pub(crate) fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

/// Atomically publishes exact bytes beneath `<artifact-root>/inputs`.
pub(crate) fn persist_artifact(
    artifact_root: &Path,
    stem: &str,
    extension: &str,
    bytes: &[u8],
) -> Result<PersistedArtifact, ArtifactError> {
    validate_fragment("stem", stem, true)?;
    validate_fragment("extension", extension, false)?;
    let digest = sha256_hex(bytes);
    let file_name = format!("{stem}-{digest}.{extension}");
    let relative_path = format!("inputs/{file_name}");
    let inputs = artifact_root.join("inputs");
    fs::create_dir_all(&inputs).map_err(|source| ArtifactError::Io {
        operation: "create artifact input directory",
        path: inputs.clone(),
        source,
    })?;
    let destination = inputs.join(file_name);

    if destination.exists() {
        verify_existing(&destination, bytes, &digest)?;
        return Ok(persisted(
            digest,
            relative_path,
            ArtifactDisposition::Reused,
        ));
    }

    let temporary = create_complete_temporary(&inputs, &digest, bytes)?;
    let disposition = match fs::hard_link(&temporary, &destination) {
        Ok(()) => ArtifactDisposition::Created,
        Err(source) if source.kind() == std::io::ErrorKind::AlreadyExists => {
            let verification = verify_existing(&destination, bytes, &digest);
            remove_published_temporary(&temporary)?;
            verification?;
            ArtifactDisposition::Reused
        }
        Err(source) => {
            remove_temporary(&temporary);
            return Err(ArtifactError::Io {
                operation: "publish artifact",
                path: destination,
                source,
            });
        }
    };
    if temporary.exists() {
        remove_published_temporary(&temporary)?;
    }
    sync_directory(&inputs)?;
    Ok(persisted(digest, relative_path, disposition))
}

/// Loads bytes only after containment and digest verification.
pub(crate) fn load_verified_artifact(
    artifact_root: &Path,
    descriptor: &ArtifactDescriptor,
) -> Result<VerifiedArtifact, ArtifactLoadError> {
    validate_descriptor(descriptor)?;
    let artifact_root =
        fs::canonicalize(artifact_root).map_err(|source| ArtifactLoadError::Io {
            operation: "resolve artifact root",
            path: artifact_root.to_path_buf(),
            source,
        })?;
    let relative = descriptor.path();
    if relative.as_os_str().is_empty()
        || relative
            .components()
            .any(|component| !matches!(component, Component::Normal(_)))
    {
        return Err(ArtifactLoadError::InvalidDescriptor {
            reason: "artifact path must be a nonempty normalized relative path".to_owned(),
        });
    }
    let unresolved = artifact_root.join(relative);
    let path = fs::canonicalize(&unresolved).map_err(|source| ArtifactLoadError::Io {
        operation: "resolve artifact",
        path: unresolved,
        source,
    })?;
    if !path.starts_with(&artifact_root) {
        return Err(ArtifactLoadError::InvalidDescriptor {
            reason: "artifact path resolves outside the artifact root".to_owned(),
        });
    }
    let bytes = fs::read(&path).map_err(|source| ArtifactLoadError::Io {
        operation: "read artifact",
        path: path.clone(),
        source,
    })?;
    let actual = sha256_hex(&bytes);
    if actual != descriptor.sha256 {
        return Err(ArtifactLoadError::DigestMismatch {
            path,
            expected: descriptor.sha256.clone(),
            actual,
        });
    }
    Ok(VerifiedArtifact { path, bytes })
}

fn persisted(sha256: String, path: String, disposition: ArtifactDisposition) -> PersistedArtifact {
    PersistedArtifact {
        descriptor: ArtifactDescriptor { sha256, path },
        disposition,
    }
}

fn validate_fragment(
    kind: &'static str,
    value: &str,
    allow_hyphen: bool,
) -> Result<(), ArtifactError> {
    let valid = !value.is_empty()
        && value.bytes().all(|byte| {
            byte.is_ascii_alphanumeric() || byte == b'_' || (allow_hyphen && byte == b'-')
        });
    if valid {
        Ok(())
    } else {
        Err(ArtifactError::InvalidFragment {
            kind,
            value: value.to_owned(),
        })
    }
}

fn validate_descriptor(descriptor: &ArtifactDescriptor) -> Result<(), ArtifactLoadError> {
    if descriptor.sha256.len() != 64
        || !descriptor
            .sha256
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(ArtifactLoadError::InvalidDescriptor {
            reason: "artifact SHA-256 must contain exactly 64 lowercase hexadecimal digits"
                .to_owned(),
        });
    }
    Ok(())
}

fn create_complete_temporary(
    directory: &Path,
    digest: &str,
    bytes: &[u8],
) -> Result<PathBuf, ArtifactError> {
    for _ in 0..1024 {
        let sequence = TEMPORARY_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let path = directory.join(format!(
            ".artifact-{digest}-{}-{sequence}.tmp",
            std::process::id()
        ));
        let mut file = match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => file,
            Err(source) if source.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(source) => {
                return Err(ArtifactError::Io {
                    operation: "create temporary artifact",
                    path,
                    source,
                });
            }
        };
        if let Err(source) = file.write_all(bytes).and_then(|()| file.sync_all()) {
            drop(file);
            remove_temporary(&path);
            return Err(ArtifactError::Io {
                operation: "write temporary artifact",
                path,
                source,
            });
        }
        return Ok(path);
    }
    Err(ArtifactError::TemporaryIdentityExhausted {
        directory: directory.to_path_buf(),
    })
}

fn verify_existing(path: &Path, expected: &[u8], digest: &str) -> Result<(), ArtifactError> {
    let actual = fs::read(path).map_err(|source| ArtifactError::Io {
        operation: "read existing artifact",
        path: path.to_path_buf(),
        source,
    })?;
    if actual == expected {
        Ok(())
    } else {
        Err(ArtifactError::DigestCollision {
            digest: digest.to_owned(),
            path: path.to_path_buf(),
        })
    }
}

fn sha256_hex(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    let mut encoded = String::with_capacity(digest.len() * 2);
    for byte in digest {
        use std::fmt::Write as _;
        write!(encoded, "{byte:02x}").expect("writing into a String cannot fail");
    }
    encoded
}

fn sync_directory(path: &Path) -> Result<(), ArtifactError> {
    File::open(path)
        .and_then(|directory| directory.sync_all())
        .map_err(|source| ArtifactError::Io {
            operation: "synchronize artifact input directory",
            path: path.to_path_buf(),
            source,
        })
}

fn remove_temporary(path: &Path) {
    let _ = fs::remove_file(path);
}

fn remove_published_temporary(path: &Path) -> Result<(), ArtifactError> {
    fs::remove_file(path).map_err(|source| ArtifactError::Io {
        operation: "remove temporary artifact",
        path: path.to_path_buf(),
        source,
    })
}

/// Failure while publishing immutable ecological input.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum ArtifactError {
    #[error("invalid artifact {kind} `{value}`")]
    InvalidFragment { kind: &'static str, value: String },
    #[error("failed to {operation} at `{path}`")]
    Io {
        operation: &'static str,
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("artifact digest collision for `{digest}` at `{path}`")]
    DigestCollision { digest: String, path: PathBuf },
    #[error("could not allocate a temporary artifact beneath `{directory}`")]
    TemporaryIdentityExhausted { directory: PathBuf },
}

/// Failure while locating or verifying immutable ecological input.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum ArtifactLoadError {
    #[error("invalid artifact descriptor: {reason}")]
    InvalidDescriptor { reason: String },
    #[error("failed to {operation} at `{path}`")]
    Io {
        operation: &'static str,
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    #[error("artifact `{path}` has SHA-256 `{actual}`, but metadata declares `{expected}`")]
    DigestMismatch {
        path: PathBuf,
        expected: String,
        actual: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn publication_is_atomic_idempotent_and_leaves_no_temporary_files() {
        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path().join("ecological-inputs");
        let bytes = b"stable ecological input";

        let created = persist_artifact(&root, "matrix", "json", bytes).unwrap();
        assert_eq!(created.disposition(), ArtifactDisposition::Created);
        let reused = persist_artifact(&root, "matrix", "json", bytes).unwrap();
        assert_eq!(reused.disposition(), ArtifactDisposition::Reused);
        assert_eq!(created.descriptor(), reused.descriptor());
        assert!(fs::read_dir(root.join("inputs")).unwrap().all(|entry| {
            !entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .ends_with(".tmp")
        }));
    }

    #[test]
    fn occupied_content_address_is_rejected_as_a_collision() {
        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path().join("ecological-inputs");
        let bytes = b"stable ecological input";
        let created = persist_artifact(&root, "matrix", "json", bytes).unwrap();
        fs::write(root.join(created.descriptor().path()), b"different").unwrap();

        assert!(matches!(
            persist_artifact(&root, "matrix", "json", bytes),
            Err(ArtifactError::DigestCollision { .. })
        ));
    }

    #[cfg(unix)]
    #[test]
    fn symlinked_artifact_cannot_escape_its_root() {
        use std::os::unix::fs::symlink;

        let temporary = tempfile::tempdir().unwrap();
        let root = temporary.path().join("ecological-inputs");
        fs::create_dir_all(root.join("inputs")).unwrap();
        let outside = temporary.path().join("outside.json");
        let bytes = b"outside";
        fs::write(&outside, bytes).unwrap();
        symlink(&outside, root.join("inputs/escaped.json")).unwrap();
        let descriptor = ArtifactDescriptor {
            sha256: sha256_hex(bytes),
            path: "inputs/escaped.json".to_owned(),
        };

        assert!(matches!(
            load_verified_artifact(&root, &descriptor),
            Err(ArtifactLoadError::InvalidDescriptor { .. })
        ));
    }
}
