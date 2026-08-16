from __future__ import annotations

from tests.qualification.capture_migration_core_operational import _capture_facts


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
