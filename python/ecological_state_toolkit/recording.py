"""Parallel, resumable conversion of Workflow recordings to NumPy arrays."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import StrEnum
from functools import partial
import hashlib
import json
import os
from pathlib import Path
import shutil
from typing import Any

import numpy as np
from scientific_workflow_reader import open_completed_recording


REQUEST_FORMAT = "ecological-state-toolkit.conversion-request.v1"
BATCH_FORMAT = "ecological-state-toolkit.processed-batch.v1"
RECORDING_FORMAT = "ecological-state-toolkit.processed-recording.v1"


class ConversionError(ValueError):
    """A request, recording, or existing processed artifact is invalid."""


class ArrayEncoding(StrEnum):
    """Supported stable ecological encodings in Workflow JSON values."""

    TENSOR_F64 = "tensor_f64"
    NONNEGATIVE_F64_VECTOR = "nonnegative_f64_vector"
    NONNEGATIVE_U32_VECTOR = "nonnegative_u32_vector"
    CATEGORICAL_LATTICE = "categorical_lattice"
    FLOAT_SCALAR = "float_scalar"
    INTEGER_SCALAR = "integer_scalar"


@dataclass(frozen=True, slots=True)
class FieldSpec:
    """One recorded field and the NPY name/encoding it should produce."""

    name: str
    encoding: ArrayEncoding
    output: str
    category_count: int | None = None

    def __post_init__(self) -> None:
        _identifier(self.name, "field name")
        _identifier(self.output, "field output")
        if self.encoding is ArrayEncoding.CATEGORICAL_LATTICE:
            if (
                isinstance(self.category_count, bool)
                or not isinstance(self.category_count, int)
                or self.category_count <= 0
            ):
                raise ValueError("categorical lattice fields require category_count > 0")
        elif self.category_count is not None:
            raise ValueError("category_count is only valid for categorical lattices")


@dataclass(frozen=True, slots=True)
class StreamSpec:
    """One Workflow stream and the fields converted from each record."""

    name: str
    fields: tuple[FieldSpec, ...]

    def __post_init__(self) -> None:
        _identifier(self.name, "stream name")
        if not self.fields:
            raise ValueError("a stream must contain at least one field")
        names = [item.name for item in self.fields]
        outputs = [item.output for item in self.fields]
        if len(names) != len(set(names)) or len(outputs) != len(set(outputs)):
            raise ValueError("stream field names and outputs must be unique")


@dataclass(frozen=True, slots=True)
class RecordingSpec:
    """One completed recording and its model-neutral conversion contract."""

    recording: Path
    identity: str
    streams: tuple[StreamSpec, ...]
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.recording, Path):
            raise TypeError("recording must be pathlib.Path")
        _nonempty_string(self.identity, "recording identity")
        if not self.streams:
            raise ValueError("a recording must contain at least one stream")
        names = [item.name for item in self.streams]
        if len(names) != len(set(names)):
            raise ValueError("recording stream names must be unique")
        try:
            json.dumps(self.metadata, sort_keys=True)
        except (TypeError, ValueError) as error:
            raise ValueError("recording metadata must be JSON-serializable") from error


@dataclass(frozen=True, slots=True)
class ArrayDescriptor:
    """A verified C-contiguous NPY array relative to the batch root."""

    path: Path
    dtype: str
    shape: tuple[int, ...]

    def as_document(self) -> dict[str, object]:
        return {"path": str(self.path), "dtype": self.dtype, "shape": list(self.shape)}


@dataclass(frozen=True, slots=True)
class ConvertedRecording:
    """Published arrays and provenance for one completed recording."""

    ordinal: int
    identity: str
    directory: Path
    source_recording: Path
    source_metadata_checksum: str
    request_checksum: str
    arrays: Mapping[str, ArrayDescriptor]
    user_metadata: Mapping[str, object]
    metadata: Mapping[str, object]

    def as_document(self) -> dict[str, object]:
        return {
            "format": RECORDING_FORMAT,
            "ordinal": self.ordinal,
            "identity": self.identity,
            "directory": str(self.directory),
            "source_recording": str(self.source_recording),
            "source_metadata_checksum": self.source_metadata_checksum,
            "request_checksum": self.request_checksum,
            "arrays": {
                name: descriptor.as_document()
                for name, descriptor in self.arrays.items()
            },
            "user_metadata": dict(self.user_metadata),
            "metadata": dict(self.metadata),
        }


ProgressCallback = Callable[[int, int, ConvertedRecording], None]


def _identifier(value: object, label: str) -> str:
    if not isinstance(value, str) or not value or any(
        character in value for character in ("/", "\\", "\0")
    ):
        raise ValueError(f"{label} must be a nonempty path-safe string")
    return value


def _nonempty_string(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a nonempty string")
    return value


def _object(value: object, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ConversionError(f"{label} must be an object")
    return dict(value)


def _array(value: object, label: str) -> list[Any]:
    if not isinstance(value, list):
        raise ConversionError(f"{label} must be an array")
    return value


def _nonnegative_integer(value: object, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ConversionError(f"{label} must be a nonnegative integer")
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return "sha256:" + digest.hexdigest()


def _checksum_document(document: Mapping[str, object]) -> str:
    encoded = json.dumps(
        document, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode()
    return "sha256:" + hashlib.sha256(encoded).hexdigest()


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        return _object(json.loads(path.read_text(encoding="utf-8")), label)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ConversionError(f"cannot read {label} at {path}") from error


def _write_json(path: Path, document: Mapping[str, object]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(document, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def _atomic_json(path: Path, document: Mapping[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    temporary.unlink(missing_ok=True)
    try:
        _write_json(temporary, document)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _decode_tensor_f64(value: object) -> np.ndarray:
    tensor = _object(value, "f64 tensor")
    if (
        tensor.get("kind") != "tensor"
        or tensor.get("version") != 1
        or tensor.get("scalar") != "f64"
    ):
        raise ConversionError("unsupported f64 tensor")
    shape = tuple(
        _nonnegative_integer(dimension, "tensor dimension")
        for dimension in _array(tensor.get("shape"), "tensor shape")
    )
    if not shape or any(dimension == 0 for dimension in shape):
        raise ConversionError("tensor shape must be nonempty and positive")
    data = np.asarray(_array(tensor.get("data"), "tensor data"), dtype=np.float64)
    if data.size != int(np.prod(shape)) or np.any(~np.isfinite(data)):
        raise ConversionError("tensor data does not match its shape")
    return np.ascontiguousarray(data.reshape(shape))


def _decode_nonnegative_f64_vector(value: object) -> np.ndarray:
    data = np.asarray(_array(value, "nonnegative f64 vector"), dtype=np.float64)
    if data.ndim != 1 or data.size == 0 or np.any(~np.isfinite(data)) or np.any(data < 0):
        raise ConversionError("value is not a finite nonnegative f64 vector")
    return np.ascontiguousarray(data)


def _decode_nonnegative_u32_vector(value: object) -> np.ndarray:
    data = np.asarray(_array(value, "nonnegative u32 vector"))
    if (
        data.ndim != 1
        or data.size == 0
        or data.dtype.kind not in "iu"
        or np.any(data < 0)
        or np.any(data > np.iinfo(np.uint32).max)
    ):
        raise ConversionError("value is not a nonnegative u32 vector")
    return np.ascontiguousarray(data, dtype=np.uint32)


def _decode_categorical_lattice(value: object, category_count: int) -> np.ndarray:
    lattice = _object(value, "categorical lattice")
    if (
        lattice.get("kind") != "square_lattice_periodic"
        or lattice.get("version") != 1
        or lattice.get("scalar") != "usize"
    ):
        raise ConversionError("unsupported categorical lattice")
    shape = tuple(
        _nonnegative_integer(dimension, "lattice dimension")
        for dimension in _array(lattice.get("shape"), "lattice shape")
    )
    if not shape or any(dimension == 0 for dimension in shape):
        raise ConversionError("lattice shape must be nonempty and positive")
    data = np.asarray(_array(lattice.get("data"), "lattice data"))
    if (
        data.size != int(np.prod(shape))
        or data.dtype.kind not in "iu"
        or np.any(data < 0)
        or np.any(data >= category_count)
    ):
        raise ConversionError("lattice values do not match shape or category count")
    maximum = category_count - 1
    dtype = np.min_scalar_type(maximum)
    if np.dtype(dtype).kind != "u":
        dtype = np.uint8
    return np.ascontiguousarray(data.reshape(shape), dtype=dtype)


def _decoder(spec: FieldSpec) -> Callable[[object], object]:
    if spec.encoding is ArrayEncoding.TENSOR_F64:
        return _decode_tensor_f64
    if spec.encoding is ArrayEncoding.NONNEGATIVE_F64_VECTOR:
        return _decode_nonnegative_f64_vector
    if spec.encoding is ArrayEncoding.NONNEGATIVE_U32_VECTOR:
        return _decode_nonnegative_u32_vector
    if spec.encoding is ArrayEncoding.CATEGORICAL_LATTICE:
        assert spec.category_count is not None
        return partial(_decode_categorical_lattice, category_count=spec.category_count)
    if spec.encoding is ArrayEncoding.FLOAT_SCALAR:
        return float
    if spec.encoding is ArrayEncoding.INTEGER_SCALAR:
        return int
    raise AssertionError(f"unhandled encoding {spec.encoding}")


def _spec_document(spec: RecordingSpec) -> dict[str, object]:
    return {
        "recording": str(spec.recording.expanduser().resolve()),
        "identity": spec.identity,
        "streams": [
            {
                "name": stream.name,
                "fields": [
                    {
                        "name": item.name,
                        "encoding": item.encoding.value,
                        "output": item.output,
                        **(
                            {"category_count": item.category_count}
                            if item.category_count is not None
                            else {}
                        ),
                    }
                    for item in stream.fields
                ],
            }
            for stream in spec.streams
        ],
        "metadata": dict(spec.metadata),
    }


def _spec_from_document(value: object) -> RecordingSpec:
    document = _object(value, "recording specification")
    streams = []
    for raw_stream in _array(document.get("streams"), "recording streams"):
        stream = _object(raw_stream, "stream specification")
        fields = []
        for raw_field in _array(stream.get("fields"), "stream fields"):
            item = _object(raw_field, "field specification")
            try:
                encoding = ArrayEncoding(item.get("encoding"))
            except ValueError as error:
                raise ConversionError("unsupported field encoding") from error
            fields.append(
                FieldSpec(
                    name=_identifier(item.get("name"), "field name"),
                    encoding=encoding,
                    output=_identifier(item.get("output"), "field output"),
                    category_count=item.get("category_count"),
                )
            )
        streams.append(
            StreamSpec(
                name=_identifier(stream.get("name"), "stream name"),
                fields=tuple(fields),
            )
        )
    recording = document.get("recording")
    if not isinstance(recording, str) or not recording:
        raise ConversionError("recording path must be a nonempty string")
    return RecordingSpec(
        recording=Path(recording),
        identity=_nonempty_string(document.get("identity"), "recording identity"),
        streams=tuple(streams),
        metadata=_object(document.get("metadata", {}), "recording metadata"),
    )


def load_request(path: Path) -> tuple[RecordingSpec, ...]:
    """Load a generic conversion request for the command-line interface."""
    if not isinstance(path, Path):
        raise TypeError("path must be pathlib.Path")
    document = _read_json(path.expanduser().resolve(), "conversion request")
    if document.get("format") != REQUEST_FORMAT:
        raise ConversionError("unsupported conversion request format")
    return tuple(
        _spec_from_document(item)
        for item in _array(document.get("recordings"), "request recordings")
    )


def _flush(array: np.ndarray) -> None:
    if isinstance(array, np.memmap):
        array.flush()


def _array_descriptor(root: Path, path: Path) -> ArrayDescriptor:
    values = np.load(path, mmap_mode="r", allow_pickle=False)
    if not values.flags.c_contiguous:
        raise ConversionError(f"processed array is not C-contiguous: {path}")
    return ArrayDescriptor(path.relative_to(root), values.dtype.str, tuple(values.shape))


def _convert_stream(reader: Any, spec: StreamSpec, directory: Path) -> None:
    count = reader.stream_record_count(spec.name)
    if count <= 0:
        raise ConversionError(f"stream {spec.name!r} is empty")
    records = iter(reader.iter_verified_records(spec.name))
    first = next(records)
    arrays: dict[str, np.memmap] = {}
    shapes: dict[str, tuple[int, ...]] = {}
    dtypes: dict[str, np.dtype[Any]] = {}
    for item in spec.fields:
        value = np.asarray(first.values[item.name])
        arrays[item.name] = np.lib.format.open_memmap(
            directory / f"{spec.name}_{item.output}.npy",
            mode="w+",
            dtype=value.dtype,
            shape=(count, *value.shape),
        )
        shapes[item.name] = value.shape
        dtypes[item.name] = value.dtype
        arrays[item.name][0] = value
    iterations = np.lib.format.open_memmap(
        directory / f"{spec.name}_iterations.npy", mode="w+", dtype=np.uint64, shape=(count,)
    )
    pending_time = directory / f".{spec.name}_physical_times.pending.npy"
    physical_times = np.lib.format.open_memmap(
        pending_time, mode="w+", dtype=np.float64, shape=(count,)
    )
    time_presence: bool | None = None

    def store_coordinate(index: int, record: Any) -> None:
        nonlocal time_presence
        iterations[index] = record.iteration
        present = record.physical_time is not None
        if time_presence is None:
            time_presence = present
        elif time_presence != present:
            raise ConversionError("physical-time presence changes within one stream")
        physical_times[index] = record.physical_time if present else np.nan

    store_coordinate(0, first)
    index = 1
    for record in records:
        if index >= count:
            raise ConversionError(f"stream {spec.name!r} has too many records")
        for item in spec.fields:
            value = np.asarray(record.values[item.name])
            if value.shape != shapes[item.name] or value.dtype != dtypes[item.name]:
                raise ConversionError(
                    f"field {item.name!r} changes shape or dtype in stream {spec.name!r}"
                )
            arrays[item.name][index] = value
        store_coordinate(index, record)
        index += 1
    if index != count:
        raise ConversionError(f"stream {spec.name!r} has too few records")
    for array in (*arrays.values(), iterations, physical_times):
        _flush(array)
    if time_presence:
        os.replace(pending_time, directory / f"{spec.name}_physical_times.npy")
    else:
        del physical_times
        pending_time.unlink()


def _document_to_result(document: Mapping[str, object]) -> ConvertedRecording:
    arrays = {
        name: ArrayDescriptor(
            path=Path(str(value["path"])),
            dtype=str(value["dtype"]),
            shape=tuple(int(item) for item in value["shape"]),
        )
        for name, value in _object(document.get("arrays"), "processed arrays").items()
    }
    return ConvertedRecording(
        ordinal=int(document["ordinal"]),
        identity=str(document["identity"]),
        directory=Path(str(document["directory"])),
        source_recording=Path(str(document["source_recording"])),
        source_metadata_checksum=str(document["source_metadata_checksum"]),
        request_checksum=str(document["request_checksum"]),
        arrays=arrays,
        user_metadata=_object(document.get("user_metadata"), "user metadata"),
        metadata=_object(document.get("metadata"), "request metadata"),
    )


def _existing_result(
    destination: Path, output_root: Path, metadata_checksum: str, request_checksum: str
) -> ConvertedRecording | None:
    path = destination / "recording.json"
    if not path.is_file():
        return None
    document = _read_json(path, "processed recording")
    if (
        document.get("format") != RECORDING_FORMAT
        or document.get("source_metadata_checksum") != metadata_checksum
        or document.get("request_checksum") != request_checksum
    ):
        return None
    for value in _object(document.get("arrays"), "processed arrays").values():
        descriptor = _object(value, "array descriptor")
        array_path = output_root / str(descriptor.get("path"))
        array = np.load(array_path, mmap_mode="r", allow_pickle=False)
        if (
            array.dtype.str != descriptor.get("dtype")
            or list(array.shape) != descriptor.get("shape")
            or not array.flags.c_contiguous
        ):
            return None
    return _document_to_result(document)


def _convert_one(job: tuple[int, RecordingSpec, Path]) -> ConvertedRecording:
    ordinal, spec, output_root = job
    recording = spec.recording.expanduser().resolve()
    metadata_path = recording / "metadata.json"
    metadata_checksum = _sha256(metadata_path)
    request_checksum = _checksum_document(_spec_document(spec))
    destination = output_root / f"member-{ordinal:06d}"
    existing = _existing_result(
        destination, output_root, metadata_checksum, request_checksum
    )
    if existing is not None:
        return existing
    if destination.exists():
        raise ConversionError(f"conflicting processed recording exists: {destination}")
    temporary = output_root / f".{destination.name}.tmp-{os.getpid()}"
    if temporary.exists():
        shutil.rmtree(temporary)
    temporary.mkdir(parents=True)
    try:
        decoders: dict[str, Callable[[object], object]] = {}
        encodings: dict[str, ArrayEncoding] = {}
        for stream in spec.streams:
            for item in stream.fields:
                previous = encodings.get(item.name)
                if previous is not None and previous is not item.encoding:
                    raise ConversionError(
                        f"field {item.name!r} has incompatible encodings across streams"
                    )
                encodings[item.name] = item.encoding
                decoders[item.name] = _decoder(item)
        reader = open_completed_recording(recording, decoders=decoders)
        for stream in spec.streams:
            _convert_stream(reader, stream, temporary)
        arrays: dict[str, ArrayDescriptor] = {}
        for path in sorted(temporary.glob("*.npy")):
            descriptor = _array_descriptor(temporary, path)
            arrays[path.stem] = ArrayDescriptor(
                Path(destination.name) / path.name,
                descriptor.dtype,
                descriptor.shape,
            )
        result = ConvertedRecording(
            ordinal=ordinal,
            identity=spec.identity,
            directory=Path(destination.name),
            source_recording=recording,
            source_metadata_checksum=metadata_checksum,
            request_checksum=request_checksum,
            arrays=arrays,
            user_metadata=dict(reader.user_metadata),
            metadata=dict(spec.metadata),
        )
        _write_json(temporary / "recording.json", result.as_document())
        os.replace(temporary, destination)
        return result
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def convert_recordings(
    recordings: Sequence[RecordingSpec],
    output_directory: Path,
    *,
    workers: int = 4,
    progress: ProgressCallback | None = None,
) -> tuple[ConvertedRecording, ...]:
    """Convert recordings concurrently with atomic, resumable member output."""
    if not isinstance(output_directory, Path):
        raise TypeError("output_directory must be pathlib.Path")
    if isinstance(workers, bool) or not isinstance(workers, int) or workers <= 0:
        raise ValueError("workers must be a positive integer")
    specs = tuple(recordings)
    if not specs:
        raise ValueError("recordings must not be empty")
    output = output_directory.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    jobs = [(ordinal, spec, output) for ordinal, spec in enumerate(specs)]
    converted: list[ConvertedRecording | None] = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs))) as pool:
        futures = {pool.submit(_convert_one, job): job[0] for job in jobs}
        completed = 0
        for future in as_completed(futures):
            ordinal = futures[future]
            result = future.result()
            converted[ordinal] = result
            completed += 1
            if progress is not None:
                progress(completed, len(jobs), result)
    if any(result is None for result in converted):
        raise ConversionError("not every recording was converted")
    return tuple(result for result in converted if result is not None)


def publish_batch_manifest(
    output_directory: Path,
    recordings: Sequence[ConvertedRecording],
) -> Path:
    """Atomically publish a generic manifest for command-line conversions."""
    if not isinstance(output_directory, Path):
        raise TypeError("output_directory must be pathlib.Path")
    output = output_directory.expanduser().resolve()
    path = output / "manifest.json"
    document = {
        "format": BATCH_FORMAT,
        "recordings": [recording.as_document() for recording in recordings],
    }
    if path.is_file():
        if _read_json(path, "batch manifest") == document:
            return path
        raise ConversionError(f"conflicting batch manifest exists: {path}")
    _atomic_json(path, document)
    return path
