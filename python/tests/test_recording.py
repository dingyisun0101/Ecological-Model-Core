from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from ecological_state_toolkit import (
    ArrayEncoding,
    FieldSpec,
    RecordingSpec,
    StreamSpec,
    convert_recordings,
)


def _recording(root: Path) -> None:
    root.mkdir()
    records = [
        {"iteration": 0, "values": [[2, 1], 3]},
        {"iteration": 1, "values": [[1, 2], 3]},
    ]
    encoded = b"".join(
        json.dumps(record, separators=(",", ":")).encode() + b"\n"
        for record in records
    )
    stream_root = root / "signal"
    stream_root.mkdir()
    chunk = stream_root / "chunk-000000.jsonl"
    chunk.write_bytes(encoded)
    metadata = {
        "format": "scientific-workflow-jsonl",
        "version": 7,
        "status": {"state": "complete"},
        "timing": {
            "created_at_utc": "2026-08-29T00:00:00Z",
            "finalized_at_utc": "2026-08-29T00:00:01Z",
            "active_duration_ns": 1,
            "continuation_count": 0,
        },
        "records": {"encoding": "json", "framing": "json_lines"},
        "time": {"iteration_name": "iteration", "physical_time_name": "physical_time"},
        "user_metadata": {"source": "test"},
        "terminal_metadata": {},
        "streams": [
            {
                "name": "signal",
                "directory": "signal",
                "sampling_interval": {"iterations": 1},
                "fields": [{"name": "abundance"}, {"name": "total"}],
                "storage": {
                    "layout": {"kind": "chunked", "target_bytes": 4096},
                    "storage_queue_bytes": 4096,
                },
                "chunks": [
                    {
                        "ordinal": 0,
                        "file": chunk.name,
                        "records": len(records),
                        "bytes": len(encoded),
                        "checksum": "sha256:" + hashlib.sha256(encoded).hexdigest(),
                        "first_iteration": 0,
                        "last_iteration": 1,
                    }
                ],
            }
        ],
    }
    (root / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")


class RecordingConversionTests(unittest.TestCase):
    def test_parallel_conversion_is_contiguous_and_resumable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            recording = root / "recording"
            _recording(recording)
            spec = RecordingSpec(
                recording=recording,
                identity="phase/task/member-a",
                streams=(
                    StreamSpec(
                        "signal",
                        (
                            FieldSpec(
                                "abundance",
                                ArrayEncoding.NONNEGATIVE_U32_VECTOR,
                                "values",
                            ),
                            FieldSpec("total", ArrayEncoding.INTEGER_SCALAR, "total"),
                        ),
                    ),
                ),
                metadata={"role": "example"},
            )
            output = root / "processed"
            first = convert_recordings((spec,), output, workers=1)
            resumed = convert_recordings((spec,), output, workers=1)
            self.assertEqual(first[0].as_document(), resumed[0].as_document())
            descriptor = first[0].arrays["signal_values"]
            values = np.load(output / descriptor.path, mmap_mode="r")
            self.assertEqual(values.dtype, np.dtype(np.uint32))
            self.assertEqual(values.shape, (2, 2))
            self.assertTrue(values.flags.c_contiguous)


if __name__ == "__main__":
    unittest.main()
