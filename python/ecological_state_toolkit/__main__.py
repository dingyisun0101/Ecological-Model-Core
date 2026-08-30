from __future__ import annotations

import argparse
from pathlib import Path

from .recording import convert_recordings, load_request, publish_batch_manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert completed ecological Workflow recordings to NPY arrays."
    )
    parser.add_argument("--request", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    arguments = parser.parse_args()

    def report(done: int, total: int, result: object) -> None:
        print(f"converted {done}/{total}: {result.identity}", flush=True)

    recordings = convert_recordings(
        load_request(arguments.request),
        arguments.output,
        workers=arguments.workers,
        progress=report,
    )
    manifest = publish_batch_manifest(arguments.output, recordings)
    print(f"published {manifest}", flush=True)


if __name__ == "__main__":
    main()
