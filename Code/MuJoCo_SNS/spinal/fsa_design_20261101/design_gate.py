"""Fail-closed all-family prerequisite for the next Hunt/curriculum campaign.

Design and individual electrical-pathway probes are allowed before release.
Existing historical runners are not silently redefined by this module. New
campaign launchers must call require_calculation_release before creating a job.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

GATE = Path(__file__).with_name("CAMPAIGN_GATE.json")
CHECKS = (
    "functional_design", "signal_units_and_domains", "time_and_body_scaling",
    "individual_pathways", "bilateral_parameterization", "control_coverage",
    "consumed_profile", "baseline_preservation",
)
OPERATIONS = {"Hunt_training", "curriculum", "whole_body_optimization"}


def _read(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _file_check(base: Path, row: dict) -> str | None:
    if not isinstance(row, dict) or not isinstance(row.get("path"), str):
        return "invalid file evidence"
    expected = row.get("sha256", "")
    if len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
        return "missing/invalid SHA256"
    path = Path(row["path"])
    if not path.is_absolute():
        path = base / path
    try:
        with path.open("rb") as stream:
            digest = hashlib.sha256()
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != expected:
            return f"changed evidence/source: {path}"
    except OSError as exc:
        return f"unreadable evidence/source: {path} ({exc.__class__.__name__})"
    return None


def inspect_gate(path: Path = GATE) -> dict:
    path = Path(path).resolve()
    gate = _read(path)
    issues = []
    if gate.get("schema") != "aarl_all_controller_fsa_design_gate_v1":
        issues.append("unknown gate schema")
    if gate.get("all_families_required_before_training") is not True:
        issues.append("all-family prerequisite is not enabled")
    families = gate.get("required_families", [])
    if not families or not all(isinstance(f, str) for f in families):
        raise ValueError("Gate must name its required controller families")
    if len(families) != len(set(families)):
        issues.append("duplicate required family")
    ledger = gate.get("family_evidence", {})
    if not isinstance(ledger, dict):
        raise ValueError("family_evidence must be an object")
    extra = sorted(set(ledger) - set(families))
    if extra:
        issues.append("unrecognized evidence families: " + ", ".join(extra))
    pending = {}
    verified = []
    for family in families:
        row = ledger.get(family)
        problems = []
        if not isinstance(row, dict):
            pending[family] = ["no reviewed family evidence"]
            continue
        if row.get("reviewed_complete") is not True:
            problems.append("family evidence has not been reviewed complete")
        checks = row.get("checks", {})
        for key in CHECKS:
            check = checks.get(key, {})
            if not isinstance(check, dict) or check.get("passed") is not True:
                problems.append(f"unverified {key}")
            elif not isinstance(check.get("evidence"), str) or not check["evidence"].strip():
                problems.append(f"no explanation/evidence for {key}")
        files = row.get("files", [])
        if not files or not isinstance(files, list):
            problems.append("no hashed design/pathway/profile/source files")
        else:
            for item in files:
                result = _file_check(path.parent, item)
                if result:
                    problems.append(result)
            kinds = {item.get("role") for item in files if isinstance(item, dict)}
            for required in ("design", "pathway_results", "applied_profile", "source"):
                if required not in kinds:
                    problems.append(f"missing hashed {required} evidence")
        if problems:
            pending[family] = problems
        else:
            verified.append(family)
    flags = ("training_allowed", "curriculum_allowed", "all_family_design_verified",
             "all_family_scaling_verified", "all_family_pathways_verified",
             "all_family_control_coverage_verified")
    closed_flags = [key for key in flags if gate.get(key) is not True]
    return {
        "machine": gate.get("machine"), "gate": str(path),
        "released": not (issues or pending or closed_flags),
        "verified_families": verified, "pending_families": pending,
        "closed_flags": closed_flags, "issues": issues,
    }


def require_calculation_release(operation: str, path: Path = GATE) -> dict:
    if operation not in OPERATIONS:
        raise ValueError(f"Unknown gated operation: {operation}")
    result = inspect_gate(path)
    if not result["released"]:
        raise RuntimeError(
            f"{operation} cannot start: all-family design/scaling/pathway/control-coverage "
            f"gate is closed ({len(result['pending_families'])} pending families; "
            f"flags={result['closed_flags']}; issues={result['issues']}). "
            f"Inspect {result['gate']}. No per-family or legacy-result bypass."
        )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", type=Path, default=GATE)
    parser.add_argument("--operation", choices=sorted(OPERATIONS))
    args = parser.parse_args()
    result = inspect_gate(args.gate)
    print(json.dumps(result, indent=2))
    return 0 if args.operation is None or result["released"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
