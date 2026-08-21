from __future__ import annotations

import json
import shutil
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier

import pytest

import chemex.baselines as baselines
from chemex.baselines import (
    ArtifactContent,
    BaselineComparisonRequest,
    BaselineLifecycleConflictError,
    BaselinePublicationIntegrityError,
    CaseDefinition,
    CaseSourceAuthority,
    CpmgBaselinePublisher,
    ExecutionSpecification,
    InputMember,
    LegacyObservationImplementation,
    Occurrence,
    ResultBundle,
    approved_scientific_anchor_case,
    cpmg_15n_ip_case,
)
from chemex.configuration.methods import Selection, read_methods
from chemex.configuration.parameters import read_defaults
from chemex.experiments.builder import build_experiments
from chemex.numerical_lanes import (
    LaneAttestation,
    LiveLaneAuthority,
    NumericalLane,
    RuntimeEnvironment,
    canonical_lanes,
)
from chemex.printers.parameters import classify_parameters, write_parameters
from chemex.runtime import AnalysisSession


def _inventory(roles: Sequence[str] = ("result:a", "result:b")) -> dict[str, object]:
    return {
        "version": "cpmg-structured-artifact-contract-v1",
        "closed": True,
        "excluded_path_components": [],
        "structured_suffixes": [".fit"],
        "required_roles": sorted(roles),
    }


def _case(*, input_hash: str = "a" * 64) -> CaseDefinition:
    return CaseDefinition.create(
        "example-anchor",
        CaseSourceAuthority("frozen-source", "b" * 64),
        {"model": "2st", "analysis": "fit"},
        (
            InputMember("parameters", input_hash, 4),
            InputMember("experiment", "c" * 64, 3),
        ),
    )


def _specification(
    case: CaseDefinition,
    *,
    workflow: Mapping[str, object] = {"kind": "fit", "method": "two-step"},
    policy: Mapping[str, object] = {"kind": "legacy-product"},
    budget: Mapping[str, object] = {"kind": "legacy-default"},
    seed: object = None,
    execution_settings: Mapping[str, object] = {"workers": "auto"},
    artifact_inventory: Mapping[str, object] | None = None,
    roles: Sequence[str] = ("Scientific anchor",),
    claims: Sequence[str] = ("legacy-observation-continuity",),
    lane_reference: str = "unqualified-local-lane-v1",
) -> ExecutionSpecification:
    return ExecutionSpecification.create(
        case,
        LegacyObservationImplementation(),
        workflow=workflow,
        lane_reference=lane_reference,
        policy=policy,
        budget=budget,
        seed=seed,
        execution_settings=execution_settings,
        artifact_inventory=artifact_inventory or _inventory(),
        roles=roles,
        claims=claims,
    )


def _live_authority(
    monkeypatch: pytest.MonkeyPatch, lane: NumericalLane
) -> LiveLaneAuthority:
    environment = RuntimeEnvironment(lane.semantics)
    monkeypatch.setattr(
        RuntimeEnvironment,
        "from_current_process",
        classmethod(lambda cls, image_digest, provenance_path=None: environment),
    )
    return lane.attest_current_process(lane.semantics.image_digest)


def _artifacts(
    roles: Sequence[str] = ("result:a", "result:b"),
) -> tuple[ArtifactContent, ...]:
    return tuple(ArtifactContent(role, role.encode("ascii")) for role in sorted(roles))


def _completed(
    case: CaseDefinition,
    specification: ExecutionSpecification,
    attempt: str,
    artifacts: Sequence[ArtifactContent] | None = None,
) -> tuple[Occurrence, ResultBundle, tuple[ArtifactContent, ...]]:
    collected = tuple(artifacts or _artifacts())
    requested = Occurrence.requested(specification, case, attempt)
    bundle = ResultBundle.create(
        requested.identity,
        specification.identity,
        specification.implementation,
        tuple(item.member for item in collected),
    )
    return requested.succeeded(bundle), bundle, collected


def _canonical_bytes(record: Mapping[str, object]) -> bytes:
    return json.dumps(
        record, ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("ascii")


def test_case_and_execution_identities_are_deterministic_and_separate() -> None:
    first = _case()
    second = CaseDefinition.create(
        "example-anchor",
        CaseSourceAuthority("frozen-source", "b" * 64),
        {"analysis": "fit", "model": "2st"},
        (
            InputMember("experiment", "c" * 64, 3),
            InputMember("parameters", "a" * 64, 4),
        ),
    )
    specification = _specification(first)

    assert first.identity == second.identity
    assert first.identity != _case(input_hash="d" * 64).identity
    assert specification.identity == _specification(first).identity
    assert specification.identity != _specification(first, seed=1).identity
    assert (
        specification.identity
        != _specification(first, budget={"kind": "other"}).identity
    )
    assert (
        specification.identity
        != _specification(first, policy={"kind": "other"}).identity
    )
    assert (
        specification.identity
        != _specification(first, execution_settings={"workers": 1}).identity
    )
    assert (
        specification.identity != _specification(first, roles=("Other role",)).identity
    )
    assert specification.identity != _specification(first, claims=("other",)).identity
    assert specification.case_identity == first.identity
    assert specification.roles == ("Scientific anchor",)
    assert (
        specification.implementation.authority_role == "LegacyObservationImplementation"
    )
    assert specification.implementation.identity != first.source_authority.identity


def test_genuine_current_process_authority_qualifies_a_matching_occurrence(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    case = _case()
    lane, other_lane = canonical_lanes()
    specification = _specification(
        case,
        lane_reference=lane.identity,
        execution_settings={"native_threads": 1, "workers": 1},
    )
    wrong_authority = _live_authority(monkeypatch, other_lane)
    authority = _live_authority(monkeypatch, lane)

    with pytest.raises(ValueError, match="requires a post-import"):
        Occurrence.requested(specification, case, "unattested")
    with pytest.raises(ValueError, match="does not match"):
        Occurrence.requested(
            specification,
            case,
            "wrong-lane",
            wrong_authority,
        )

    occurrence = Occurrence.requested(specification, case, "qualified", authority)
    evidence = LaneAttestation.from_record(authority.to_record())

    assert occurrence.lane_reference == lane.identity
    assert occurrence.lane_attestation_identity == evidence.identity
    assert Occurrence.from_record(occurrence.to_record(), specification) == occurrence
    with pytest.raises(ValueError, match="cannot carry lane authority"):
        Occurrence.requested(_specification(case), case, "retroactive", authority)
    incompatible_specification = _specification(
        case,
        lane_reference=lane.identity,
        execution_settings={"native_threads": 1, "workers": 8},
    )
    with pytest.raises(ValueError, match="lane concurrency"):
        Occurrence.requested(
            incompatible_specification, case, "wrong-workers", authority
        )

    unqualified_specification = _specification(case)
    unqualified = Occurrence.requested(
        unqualified_specification, case, "legacy-wire-format"
    )
    legacy_record = unqualified.to_record()
    assert "lane_attestation_identity" not in legacy_record
    assert (
        Occurrence.from_record(legacy_record, unqualified_specification) == unqualified
    )
    legacy_record["lane_attestation_identity"] = None
    with pytest.raises(ValueError, match="unknown or missing fields"):
        Occurrence.from_record(legacy_record, unqualified_specification)


def test_deserialized_evidence_cannot_qualify_a_new_occurrence() -> None:
    case = _case()
    lane_identity = "d" * 64
    specification = _specification(
        case,
        lane_reference=lane_identity,
        execution_settings={"native_threads": 1, "workers": 1},
    )
    fabricated = LaneAttestation(
        lane_identity,
        "e" * 64,
        1,
        1,
        "POST_IMPORT_CURRENT_PROCESS",
    )
    evidence = LaneAttestation.from_record(fabricated.to_record())

    with pytest.raises(TypeError, match="live current-process lane authority"):
        Occurrence.requested(specification, case, "fabricated", evidence)


@pytest.mark.parametrize(
    ("field", "values"),
    (
        ("roles", ("Scientific anchor", "A")),
        ("roles", ("Scientific anchor", "Scientific anchor")),
        ("claims", ("z", "a")),
        ("claims", ("a", "a")),
    ),
)
def test_execution_specification_create_rejects_noncanonical_roles_and_claims(
    field: str, values: Sequence[str]
) -> None:
    with pytest.raises(ValueError, match="canonically"):
        if field == "roles":
            _specification(_case(), roles=values)
        else:
            _specification(_case(), claims=values)


def test_records_reject_noncanonical_sequences_and_duplicate_members() -> None:
    case = _case()
    specification = _specification(
        case, roles=("A", "Scientific anchor"), claims=("a", "z")
    )
    occurrence, bundle, _ = _completed(case, specification, "attempt-a")

    reversed_case = case.to_record()
    raw_inputs = reversed_case["inputs"]
    assert isinstance(raw_inputs, list)
    reversed_case["inputs"] = list(reversed(raw_inputs))
    with pytest.raises(ValueError, match="canonically"):
        CaseDefinition.from_record(reversed_case)

    reversed_bundle = bundle.to_record()
    raw_members = reversed_bundle["members"]
    assert isinstance(raw_members, list)
    reversed_bundle["members"] = list(reversed(raw_members))
    with pytest.raises(ValueError, match="canonically"):
        ResultBundle.from_record(reversed_bundle, occurrence)

    duplicate_bundle = bundle.to_record()
    duplicate_members = duplicate_bundle["members"]
    assert isinstance(duplicate_members, list)
    duplicate_bundle["members"] = [duplicate_members[0], duplicate_members[0]]
    with pytest.raises(ValueError, match="unique"):
        ResultBundle.from_record(duplicate_bundle, occurrence)

    reordered_specification = specification.to_record()
    reordered_specification["roles"] = ["Scientific anchor", "A"]
    with pytest.raises(ValueError, match="canonically"):
        ExecutionSpecification.from_record(reordered_specification)


def test_result_bundle_record_requires_the_exact_successful_occurrence() -> None:
    case = _case()
    specification = _specification(case)
    successful, bundle, artifacts = _completed(case, specification, "attempt-a")
    requested = Occurrence.requested(specification, case, "attempt-a")
    failed = requested.failed("RunnerError")
    other_bundle = ResultBundle.create(
        requested.identity,
        specification.identity,
        specification.implementation,
        tuple(
            ArtifactContent(artifact.role, b"other").member for artifact in artifacts
        ),
    )

    with pytest.raises(ValueError, match="successful occurrence"):
        ResultBundle.from_record(bundle.to_record(), requested)
    with pytest.raises(ValueError, match="successful occurrence"):
        ResultBundle.from_record(bundle.to_record(), failed)
    with pytest.raises(ValueError, match="does not belong"):
        ResultBundle.from_record(bundle.to_record(), requested.succeeded(other_bundle))
    assert ResultBundle.from_record(bundle.to_record(), successful) == bundle


def test_changed_artifact_bytes_change_member_and_bundle_identities() -> None:
    case = _case()
    specification = _specification(case)
    requested = Occurrence.requested(specification, case, "attempt-a")
    original = ArtifactContent("result:a", b"original")
    changed = ArtifactContent("result:a", b"changed")
    shared = ArtifactContent("result:b", b"shared")
    first = ResultBundle.create(
        requested.identity,
        specification.identity,
        specification.implementation,
        (original.member, shared.member),
    )
    second = ResultBundle.create(
        requested.identity,
        specification.identity,
        specification.implementation,
        (changed.member, shared.member),
    )

    assert original.member.content_hash != changed.member.content_hash
    assert first.manifest_identity != second.manifest_identity
    assert first.identity != second.identity
    reordered_specification = specification.to_record()
    reordered_specification["claims"] = ["z", "a"]
    with pytest.raises(ValueError, match="canonically"):
        ExecutionSpecification.from_record(reordered_specification)


def test_occurrence_reservation_allows_one_terminal_result_and_distinct_reruns(
    tmp_path: Path,
) -> None:
    case = _case()
    specification = _specification(case)
    publisher = CpmgBaselinePublisher(tmp_path / "evidence")
    successful, bundle, artifacts = _completed(case, specification, "attempt-a")
    requested = Occurrence.requested(specification, case, "attempt-a")
    publisher.reserve(case, specification, requested)
    publisher.publish(case, specification, successful, bundle, artifacts)

    changed_artifacts = (
        ArtifactContent("result:a", b"changed"),
        ArtifactContent("result:b", b"result:b"),
    )
    changed_bundle = ResultBundle.create(
        requested.identity,
        specification.identity,
        specification.implementation,
        tuple(item.member for item in changed_artifacts),
    )
    with pytest.raises(BaselineLifecycleConflictError) as error:
        publisher.publish(
            case,
            specification,
            requested.succeeded(changed_bundle),
            changed_bundle,
            changed_artifacts,
        )
    assert error.value.occurrence == successful
    assert publisher.terminal_occurrence(case, specification, requested) == successful

    failed_request = Occurrence.requested(specification, case, "attempt-b")
    publisher.reserve(case, specification, failed_request)
    publisher.record_failure(case, specification, failed_request, "RunnerError")
    failed_bundle = ResultBundle.create(
        failed_request.identity,
        specification.identity,
        specification.implementation,
        tuple(item.member for item in artifacts),
    )
    with pytest.raises(BaselineLifecycleConflictError) as error:
        publisher.publish(
            case,
            specification,
            failed_request.succeeded(failed_bundle),
            failed_bundle,
            artifacts,
        )
    assert error.value.occurrence == publisher.terminal_occurrence(
        case, specification, failed_request
    )

    rerun = Occurrence.requested(specification, case, "attempt-c")
    publisher.reserve(case, specification, rerun)
    assert rerun.identity != requested.identity
    assert rerun.case_identity == requested.case_identity
    assert (
        rerun.execution_specification_identity
        == requested.execution_specification_identity
    )


def test_terminal_conflicts_expose_only_the_persisted_occurrence(
    tmp_path: Path,
) -> None:
    case = _case()
    specification = _specification(case)
    publisher = CpmgBaselinePublisher(tmp_path / "evidence")

    success, bundle, artifacts = _completed(case, specification, "success-wins")
    requested = Occurrence.requested(specification, case, "success-wins")
    publisher.reserve(case, specification, requested)
    publisher.publish(case, specification, success, bundle, artifacts)
    with pytest.raises(BaselineLifecycleConflictError) as error:
        publisher.record_failure(case, specification, requested, "loser")
    assert error.value.occurrence == success
    assert publisher.terminal_occurrence(case, specification, requested) == success

    requested = Occurrence.requested(specification, case, "failure-wins")
    publisher.reserve(case, specification, requested)
    failed = publisher.record_failure(case, specification, requested, "winner")
    bundle = ResultBundle.create(
        requested.identity,
        specification.identity,
        specification.implementation,
        tuple(item.member for item in artifacts),
    )
    with pytest.raises(BaselineLifecycleConflictError) as error:
        publisher.publish(
            case, specification, requested.succeeded(bundle), bundle, artifacts
        )
    assert error.value.occurrence == failed


def test_concurrent_success_and_failure_expose_one_persisted_winner(
    tmp_path: Path,
) -> None:
    case = _case()
    specification = _specification(case)
    publisher = CpmgBaselinePublisher(tmp_path / "evidence")
    requested = Occurrence.requested(specification, case, "success-versus-failure")
    publisher.reserve(case, specification, requested)
    successful, bundle, artifacts = _completed(
        case, specification, "success-versus-failure"
    )
    barrier = Barrier(2)

    def publish_success() -> baselines.PublishedEvidence:
        barrier.wait()
        return publisher.publish(case, specification, successful, bundle, artifacts)

    def publish_failure() -> Occurrence:
        barrier.wait()
        return publisher.record_failure(case, specification, requested, "loser")

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(publish_success), executor.submit(publish_failure)]
    results = [future.exception() or future.result() for future in futures]
    persisted = publisher.terminal_occurrence(case, specification, requested)
    winners = [result for result in results if not isinstance(result, Exception)]
    conflicts = [
        result
        for result in results
        if isinstance(result, BaselineLifecycleConflictError)
    ]

    assert len(winners) == 1
    assert len(conflicts) == 1
    assert conflicts[0].occurrence == persisted
    if isinstance(winners[0], baselines.PublishedEvidence):
        assert persisted == successful
        assert publisher.read(bundle.identity).occurrence == persisted
    else:
        assert winners[0] == persisted
        assert persisted.lifecycle == "FAILED"
        assert not tuple(publisher.root.glob("attempts/*/terminal/success"))


def test_invalid_terminal_collision_is_an_integrity_failure(tmp_path: Path) -> None:
    case = _case()
    specification = _specification(case)
    publisher = CpmgBaselinePublisher(tmp_path / "evidence")
    requested = Occurrence.requested(specification, case, "malformed-terminal")
    publisher.reserve(case, specification, requested)
    terminal = publisher._attempt_directory(requested) / "terminal"
    terminal.mkdir()
    (terminal / "malformed").write_text("not a terminal record", encoding="ascii")
    successful, bundle, artifacts = _completed(
        case, specification, "malformed-terminal"
    )

    with pytest.raises(BaselinePublicationIntegrityError):
        publisher.publish(case, specification, successful, bundle, artifacts)


def test_concurrent_success_contenders_install_one_terminal_bundle(
    tmp_path: Path,
) -> None:
    case = _case()
    specification = _specification(case)
    publisher = CpmgBaselinePublisher(tmp_path / "evidence")
    requested = Occurrence.requested(specification, case, "concurrent")
    publisher.reserve(case, specification, requested)
    barrier = Barrier(4)

    def contend(index: int) -> object:
        artifacts = (
            ArtifactContent("result:a", f"a-{index}".encode()),
            ArtifactContent("result:b", f"b-{index}".encode()),
        )
        bundle = ResultBundle.create(
            requested.identity,
            specification.identity,
            specification.implementation,
            tuple(item.member for item in artifacts),
        )
        barrier.wait()
        return publisher.publish(
            case, specification, requested.succeeded(bundle), bundle, artifacts
        )

    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(contend, index) for index in range(4)]
    results = [future.exception() or future.result() for future in futures]
    winners = [result for result in results if not isinstance(result, Exception)]
    assert len(winners) == 1
    assert all(
        isinstance(result, BaselineLifecycleConflictError)
        for result in results
        if isinstance(result, Exception)
    )
    winner = winners[0]
    assert isinstance(winner, baselines.PublishedEvidence)
    persisted = publisher.terminal_occurrence(case, specification, requested)
    assert persisted == winner.occurrence
    assert all(
        result.occurrence == persisted
        for result in results
        if isinstance(result, BaselineLifecycleConflictError)
    )


def test_publication_manifest_is_canonical_and_identifies_the_terminal_state(
    tmp_path: Path,
) -> None:
    case = _case()
    specification = _specification(case)
    occurrence, bundle, artifacts = _completed(case, specification, "attempt-a")
    publisher = CpmgBaselinePublisher(tmp_path / "evidence")
    publisher.reserve(
        case, specification, Occurrence.requested(specification, case, "attempt-a")
    )
    published = publisher.publish(case, specification, occurrence, bundle, artifacts)
    assert publisher.read(bundle.identity) == published

    manifest_path = published.location / "manifest.json"
    canonical = manifest_path.read_bytes()
    manifest = json.loads(canonical)
    assert manifest["manifest_identity"] == published.manifest_identity

    manifest_path.write_bytes(json.dumps(manifest, indent=2).encode("ascii"))
    with pytest.raises(ValueError, match="canonically"):
        publisher.read(bundle.identity)

    for field, replacement in (
        ("member_manifest_identity", "0" * 64),
        ("manifest_identity", "0" * 64),
    ):
        manifest_path.write_bytes(canonical)
        tampered = json.loads(canonical)
        tampered[field] = replacement
        manifest_path.write_bytes(_canonical_bytes(tampered))
        with pytest.raises(ValueError, match="identity"):
            publisher.read(bundle.identity)

    manifest_path.write_bytes(canonical)
    tampered = json.loads(canonical)
    raw_occurrence = tampered["occurrence"]
    assert isinstance(raw_occurrence, dict)
    raw_occurrence["lifecycle"] = "FAILED"
    raw_occurrence["result_bundle_identity"] = None
    raw_occurrence["failure_code"] = "RunnerError"
    manifest_path.write_bytes(_canonical_bytes(tampered))
    with pytest.raises(ValueError):
        publisher.read(bundle.identity)

    manifest_path.write_bytes(canonical)
    tampered = json.loads(canonical)
    raw_bundle = tampered["bundle"]
    assert isinstance(raw_bundle, dict)
    raw_bundle["identity"] = "0" * 64
    manifest_path.write_bytes(_canonical_bytes(tampered))
    with pytest.raises(ValueError, match="missing"):
        publisher.read(bundle.identity)

    manifest_path.write_bytes(canonical[:-1] + b',"case":null}')
    with pytest.raises(ValueError):
        publisher.read(bundle.identity)

    manifest_path.write_bytes(canonical)
    member_path = published.location / "members" / bundle.members[0].content_hash
    member_path.write_bytes(b"tampered member")
    with pytest.raises(ValueError, match="hash mismatch"):
        publisher.read(bundle.identity)


def test_comparison_request_keeps_comparison_identity_out_of_anchor_execution() -> None:
    case = _case()
    specification = _specification(case)
    first_occurrence, first, _ = _completed(case, specification, "attempt-a")
    _second_occurrence, second, _ = _completed(case, specification, "attempt-b")
    request = BaselineComparisonRequest.create(
        first,
        second,
        "NUMERICAL_ARTIFACT_COMPARISON",
        {"version": "comparison-policy-v1"},
    )

    assert first_occurrence.lifecycle == "SUCCEEDED"
    assert BaselineComparisonRequest.from_record(request.to_record()) == request


def _cpmg_specification(anchor: Path) -> tuple[CaseDefinition, ExecutionSpecification]:
    inputs = baselines._capture_cpmg_inputs(anchor)
    case = baselines._case_from_cpmg_inputs(inputs)
    return case, baselines.cpmg_15n_ip_legacy_specification(
        case,
        artifact_inventory=baselines._cpmg_artifact_inventory(inputs),
        implementation=LegacyObservationImplementation(),
    )


def _write_complete_cpmg_output(
    output: Path, specification: ExecutionSpecification
) -> None:
    for role in baselines.CpmgBaselinePublisher._required_roles(specification):
        relative = role.removeprefix("legacy-output:")
        path = output / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(role.encode("ascii"))


def _write_complete_approved_anchor_output(
    output: Path, snapshot: Path, anchor_name: baselines.ApprovedAnchorName
) -> None:
    anchor = baselines._approved_anchor(anchor_name)
    inputs = baselines._capture_anchor_inputs(anchor, snapshot)
    inventory = baselines._anchor_artifact_inventory(anchor, inputs)
    roles = inventory["required_roles"]
    assert isinstance(roles, list)
    for role in roles:
        assert isinstance(role, str)
        if not role.startswith("legacy-output:"):
            continue
        relative = role.removeprefix("legacy-output:")
        path = output / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(role.encode("ascii"))


def test_cpmg_closed_artifact_contract_rejects_partial_and_unexpected_output(
    tmp_path: Path,
) -> None:
    _case, specification = _cpmg_specification(Path("examples/Experiments/CPMG_15N_IP"))
    output = tmp_path / "Output"
    output.mkdir()
    (output / "one.fit").write_bytes(b"one")
    with pytest.raises(ValueError, match="closed required"):
        baselines._legacy_result_artifacts(output, specification)

    shutil.rmtree(output)
    output.mkdir()
    _write_complete_cpmg_output(output, specification)
    artifacts = baselines._legacy_result_artifacts(output, specification)
    assert len(artifacts) == 340

    missing = output / artifacts[0].role.removeprefix("legacy-output:")
    missing.unlink()
    with pytest.raises(ValueError, match="closed required"):
        baselines._legacy_result_artifacts(output, specification)

    missing.write_bytes(b"restored")
    (output / "unexpected.toml").write_bytes(b"unexpected")
    with pytest.raises(ValueError, match="closed required"):
        baselines._legacy_result_artifacts(output, specification)


@pytest.mark.parametrize(
    ("anchor_name", "directory"),
    (
        ("cest-13c-label-cn", "examples/Experiments/CEST_13C_LABEL_CN"),
        ("2st-binding", "examples/Combinations/2stBinding"),
        ("dcest-fifu-drd", "examples/Experiments/DCEST_15N_3States"),
    ),
)
def test_new_anchor_artifact_contracts_are_closed(
    tmp_path: Path,
    anchor_name: baselines.ApprovedAnchorName,
    directory: str,
) -> None:
    anchor_directory = Path(directory)
    anchor = baselines._approved_anchor(anchor_name)
    inputs = baselines._capture_anchor_inputs(anchor, anchor_directory)
    case = baselines._case_from_anchor_inputs(anchor, inputs)
    inventory = baselines._anchor_artifact_inventory(anchor, inputs)
    specification = baselines._anchor_legacy_specification(
        anchor,
        case,
        artifact_inventory=inventory,
        implementation=LegacyObservationImplementation(),
    )
    output = tmp_path / anchor.output_directory
    _write_complete_approved_anchor_output(output, anchor_directory, anchor_name)

    artifacts = baselines._legacy_result_artifacts(output, specification)

    assert len(artifacts) > 1
    missing = output / artifacts[0].role.removeprefix("legacy-output:")
    missing.unlink()
    with pytest.raises(ValueError, match="closed required"):
        baselines._legacy_result_artifacts(output, specification)

    missing.write_bytes(b"restored")
    (output / "unexpected.toml").write_bytes(b"unexpected")
    with pytest.raises(ValueError, match="closed required"):
        baselines._legacy_result_artifacts(output, specification)


def test_cest_anchor_omits_only_the_empty_step1_fixed_report(tmp_path: Path) -> None:
    anchor_directory = Path("examples/Experiments/CEST_13C_LABEL_CN")
    anchor = baselines._approved_anchor("cest-13c-label-cn")
    inputs = baselines._capture_anchor_inputs(anchor, anchor_directory)
    case = baselines._case_from_anchor_inputs(anchor, inputs)
    inventory = baselines._anchor_artifact_inventory(anchor, inputs)
    specification = baselines._anchor_legacy_specification(
        anchor,
        case,
        artifact_inventory=inventory,
        implementation=LegacyObservationImplementation(),
    )
    required_roles = inventory["required_roles"]
    assert isinstance(required_roles, list)

    session = AnalysisSession.create()
    session.set_model(anchor.model)
    experiments = build_experiments(
        sorted(anchor_directory.glob("Experiments/*.toml")),
        Selection(include=None, exclude=None),
        session=session,
    )
    session.parameters.set_defaults(
        read_defaults([anchor_directory / "Parameters/parameters.toml"])
    )
    step1 = read_methods([anchor_directory / "Methods/method.toml"])["STEP1"]
    experiments.select(step1.selection)
    session.parameters.set_parameter_status(step1)
    classified = classify_parameters(experiments, session.parameters)
    assert not classified.fixed

    parameter_output = tmp_path / "classified"
    write_parameters(experiments, parameter_output, session.parameters)
    assert not (parameter_output / "Parameters/fixed.toml").exists()
    assert "legacy-output:STEP1/Parameters/fixed.toml" not in required_roles

    output = tmp_path / anchor.output_directory
    _write_complete_approved_anchor_output(
        output, anchor_directory, "cest-13c-label-cn"
    )
    artifacts = baselines._legacy_result_artifacts(output, specification)
    assert len(artifacts) == 90

    unexpected = output / "STEP1/Parameters/fixed.toml"
    unexpected.write_bytes(b"")
    with pytest.raises(ValueError, match="closed required"):
        baselines._legacy_result_artifacts(output, specification)


def test_dcest_anchor_derives_decimal_output_names_from_legacy_suffix_rules() -> None:
    anchor_directory = Path("examples/Experiments/DCEST_15N_3States")
    anchor = baselines._approved_anchor("dcest-fifu-drd")
    inputs = baselines._capture_anchor_inputs(anchor, anchor_directory)
    inventory = baselines._anchor_artifact_inventory(anchor, inputs)
    required_roles = inventory["required_roles"]
    assert isinstance(required_roles, list)

    decimal_inputs = {
        item.relative_path
        for item in anchor.experiments
        if item.relative_path in {"Experiments/1.25hz.toml", "Experiments/7.5hz.toml"}
    }
    assert decimal_inputs == {
        "Experiments/1.25hz.toml",
        "Experiments/7.5hz.toml",
    }
    expected = {
        f"legacy-output:Data/{Path(path).with_suffix('').with_suffix('.dat').name}"
        for path in decimal_inputs
    } | {
        f"legacy-output:Plots/{Path(path).with_suffix('').with_suffix('.fit').name}"
        for path in decimal_inputs
    }
    assert expected == {
        "legacy-output:Data/1.dat",
        "legacy-output:Data/7.dat",
        "legacy-output:Plots/1.fit",
        "legacy-output:Plots/7.fit",
    }
    assert expected <= set(required_roles)
    assert not {
        "legacy-output:Data/1.25hz.dat",
        "legacy-output:Data/7.5hz.dat",
        "legacy-output:Plots/1.25hz.fit",
        "legacy-output:Plots/7.5hz.fit",
    } & set(required_roles)


@pytest.mark.parametrize(
    ("anchor_name", "directory", "case_name"),
    (
        ("cpmg-15n-ip", "examples/Experiments/CPMG_15N_IP", "cpmg-15n-ip"),
        (
            "cest-13c-label-cn",
            "examples/Experiments/CEST_13C_LABEL_CN",
            "cest-13c-label-cn",
        ),
        ("2st-binding", "examples/Combinations/2stBinding", "2st-binding"),
        (
            "dcest-fifu-drd",
            "examples/Experiments/DCEST_15N_3States",
            "dcest-fifu-drd",
        ),
    ),
)
def test_exactly_the_four_approved_shipped_anchor_cases_are_frozen(
    anchor_name: baselines.ApprovedAnchorName, directory: str, case_name: str
) -> None:
    case = approved_scientific_anchor_case(anchor_name, Path(directory))

    assert case.name == case_name
    assert case.source_authority == CaseSourceAuthority(
        "d5ed0c87e8ce7a7f17745feea346af4dfbae7ecf",
        "f0fb2ffc7b1a5ecd1bf7ac43956fc4861b96c058d158948b68b4e97027a6086a",
    )
    if anchor_name == "cpmg-15n-ip":
        assert case.identity == (
            "18d9a4fb78e300a474132313c898c78240be7aa6ae7eb2f3385d2fed69e6ebd3"
        )
        assert case == cpmg_15n_ip_case(Path(directory))


def test_unapproved_shipped_anchor_name_is_rejected() -> None:
    with pytest.raises(ValueError, match="approved scientific anchor"):
        approved_scientific_anchor_case(
            "cest-15n-label-cn", Path("examples/Experiments/CEST_15N_LABEL_CN")
        )
