from __future__ import annotations

from pathlib import Path

from chemex.migration_core_operational import (
    OPERATIONAL_REQUIREMENTS,
    OperationalReplayCapture,
    eligible_operational_requirements,
)
from tests.qualification.capture_migration_core_operational import _capture_facts

_REPOSITORY_ROOT = Path(__file__).parents[1]


def test_operational_capture_observes_required_runtime_facts() -> None:
    facts = _capture_facts()

    assert all(facts["serialization"].values())
    assert all(facts["multiprocessing"].values())
    assert facts["cache"] == {
        "hits": 1,
        "misses": 2,
        "changed_frame_invalidated": True,
    }
    assert all(facts["deterministic_replay"].values())
    assert all(facts["fail_closed"].values())


def test_canonical_operational_capture_reconstructs_exact_qualified_scope() -> None:
    capture = OperationalReplayCapture.from_bytes(
        (
            _REPOSITORY_ROOT
            / "tests/fixtures/migration_core_operational_replay_v1.json"
        ).read_bytes()
    )

    assert capture.identity == (
        "8a0394e50bef697dce880ee35bf7d77d2f8ffa3c8fc8e81ab0a6afbb4d15a59a"
    )
    assert capture.occurrence.lifecycle == "SUCCEEDED"
    assert (
        eligible_operational_requirements(
            capture,
            source_commit="fb156f86f431a90c65d1a7285bdb6532ab2c51ec",
            lockfile_hash=(
                "cc7a8e08d8fb8f1ea4255b63452598f6dbe041a8b4024de0f3af065020088004"
            ),
            lane_identity=(
                "4f1e7dab3384f88149fd33befb09ba3e96730dc336e9427b224a39a8a7167e4f"
            ),
            attestation_identity=(
                "3b3e0bc184826d61ec6652194486c907a5faa4c64b68ec03bbe60b63c660d687"
            ),
            environment_identity=(
                "cc5359f90df35ec9b60fd56e483911745209519ecce37d69e17c4edd6ea3604f"
            ),
        )
        == OPERATIONAL_REQUIREMENTS
    )
