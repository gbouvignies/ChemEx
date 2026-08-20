"""Diagnostic #652 proof that a separately identified wheel runs in one lane."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import cast

import chemex.messages as application_module
from chemex.baselines import LegacyObservationImplementation
from chemex.numerical_lanes import LaneRole, NumericalLane, RuntimeEnvironment


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", required=True)
    parser.add_argument(
        "--role",
        required=True,
        choices=("CANONICAL_NUMERICAL", "PYTHON_COMPATIBILITY"),
    )
    parser.add_argument("--image-digest", required=True)
    return parser


def main() -> None:
    args = _parser().parse_args()
    lane = NumericalLane.capture_current_process(
        args.name,
        cast("LaneRole", args.role),
        args.image_digest,
    )
    environment = RuntimeEnvironment.from_current_process(args.image_digest)
    authority = lane.attest_current_process(args.image_digest)
    implementation = LegacyObservationImplementation.from_current_package()
    module_path = Path(application_module.__file__).resolve()
    print(
        json.dumps(
            {
                "application_module_sha256": hashlib.sha256(
                    module_path.read_bytes()
                ).hexdigest(),
                "application_module_path": str(module_path),
                "attestation": authority.to_record(),
                "environment": environment.to_record(),
                "implementation": implementation.to_record(),
                "implementation_wheel_sha256": os.environ[
                    "CHEMEX_IMPLEMENTATION_WHEEL_SHA256"
                ],
                "lane": lane.to_record(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
