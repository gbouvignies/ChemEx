"""Render current covariance and constrained-uncertainty evidence."""

from __future__ import annotations

import json
from pathlib import Path
from typing import cast

from chemex.optimize.uncertainty import (
    RootAnchoredBlockCovarianceEvidence,
    UncertaintyEvidence,
)
from chemex.printers.parameters import uncertainty_unavailable_reason


def write_json(path: Path, record: object) -> None:
    """Write one deterministic, strict JSON evidence artifact."""
    path.write_text(
        json.dumps(
            record,
            allow_nan=False,
            ensure_ascii=True,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def write_uncertainty(path: Path, uncertainty: UncertaintyEvidence | None) -> None:
    """Write covariance and constrained evidence in distinct typed families."""
    if uncertainty is None:
        return
    record = uncertainty.to_record()
    payload = record.get("payload")
    if not isinstance(payload, dict):
        raise TypeError("Uncertainty evidence has no typed payload")
    common = {
        "schema_version": record["schema_version"],
        "accepted_result_identity": uncertainty.accepted_result_identity,
        "accepted_occurrence_identity": uncertainty.accepted_occurrence_identity,
        "request_identity": uncertainty.request_identity,
        "policy_identity": uncertainty.policy_identity,
        "bundle_identity": uncertainty.identity,
    }
    covariance_path = path / "Covariance"
    covariance_path.mkdir()
    write_json(
        covariance_path / "evidence.json",
        common
        | {
            "artifact_type": "native_covariance_evidence",
            "residual_jacobian": payload.get("residual_jacobian"),
            "rank_diagnostic": payload.get("rank_diagnostic"),
            "covariance": payload.get("covariance"),
            "marginal_standard_errors": payload.get("marginal_errors"),
            "correlations": payload.get("correlations"),
            "operations": payload.get("operations"),
            "failures": payload.get("failures"),
        },
    )
    if uncertainty.requested_output_scope:
        constrained_path = path / "Constrained"
        constrained_path.mkdir()
        write_json(
            constrained_path / "evidence.json",
            common
            | {
                "artifact_type": "native_constrained_evidence",
                "requested_output_scope": list(uncertainty.requested_output_scope),
                "constraint_jacobian": payload.get("constraint_jacobian"),
                "constrained_propagation": payload.get("constrained_propagation"),
                "marginal_standard_errors": payload.get("constrained_marginal_errors"),
                "correlations": payload.get("constrained_correlations"),
                "operations": payload.get("operations"),
                "failures": payload.get("failures"),
            },
        )


def write_block_uncertainty(
    path: Path,
    evidence: RootAnchoredBlockCovarianceEvidence | None,
) -> None:
    """Write root-authoritative failure-isolated covariance blocks."""
    if evidence is None:
        return
    covariance_path = path / "Covariance"
    covariance_path.mkdir(exist_ok=True)
    record = evidence.to_record()
    blocks = record["blocks"]
    if not isinstance(blocks, list):
        raise TypeError("Block covariance record has an invalid payload")
    for block, block_record in zip(evidence.blocks, blocks, strict=True):
        if not isinstance(block_record, dict):
            raise TypeError("Block covariance record has an invalid block")
        typed_block_record = cast("dict[str, object]", block_record)
        typed_block_record["unavailable_reason"] = (
            ""
            if block.unavailable_kind is None
            else uncertainty_unavailable_reason(block.unavailable_kind)
        )
    write_json(covariance_path / "blocks.json", record)
