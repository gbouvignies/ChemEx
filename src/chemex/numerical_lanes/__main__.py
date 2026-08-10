"""Build-time capture and post-import attestation commands for #588 lanes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import cast

from . import LaneRole, NumericalLane, _strict_json_loads


def _lane(path: Path) -> NumericalLane:
    raw = _strict_json_loads(path.read_text(encoding="ascii"))
    if not isinstance(raw, dict):
        raise TypeError("Numerical-lane manifest must be a record")
    return NumericalLane.from_record(cast("dict[str, object]", raw))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    capture = subparsers.add_parser("capture")
    capture.add_argument("--name", required=True)
    capture.add_argument(
        "--role",
        required=True,
        choices=("CANONICAL_NUMERICAL", "PYTHON_COMPATIBILITY"),
    )
    capture.add_argument("--image-digest", required=True)

    attest = subparsers.add_parser("attest")
    attest.add_argument("--manifest", type=Path, required=True)
    attest.add_argument("--image-digest", required=True)

    compatibility = subparsers.add_parser("validate-compatibility")
    compatibility.add_argument("--canonical", type=Path, required=True)
    compatibility.add_argument("--compatibility", type=Path, required=True)
    return parser


def main() -> None:
    args = _parser().parse_args()
    if args.command == "capture":
        lane = NumericalLane.capture_current_process(
            args.name,
            cast("LaneRole", args.role),
            args.image_digest,
        )
        print(json.dumps(lane.to_record(), indent=2, sort_keys=True))
    elif args.command == "attest":
        lane = _lane(args.manifest)
        attestation = lane.attest_current_process(args.image_digest)
        print(json.dumps(attestation.to_record(), indent=2, sort_keys=True))
    else:
        canonical = _lane(args.canonical)
        compatibility = _lane(args.compatibility)
        canonical.validate_compatibility_lane(compatibility)
        print(json.dumps(canonical.compatibility_delta(compatibility), sort_keys=True))


if __name__ == "__main__":
    main()
