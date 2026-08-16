"""Attest the canonical lane before loading the operational capture workload."""

from __future__ import annotations

import argparse
from pathlib import Path

from chemex.numerical_lanes import canonical_lanes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--lockfile-hash", required=True)
    parser.add_argument("--image-digest", required=True)
    arguments = parser.parse_args()
    authority = canonical_lanes()[0].attest_current_process(arguments.image_digest)

    from tests.qualification.capture_migration_core_operational import capture

    arguments.output.write_bytes(
        capture(
            source_commit=arguments.source_commit,
            lockfile_hash=arguments.lockfile_hash,
            authority=authority,
        ).to_bytes()
    )


if __name__ == "__main__":
    main()
