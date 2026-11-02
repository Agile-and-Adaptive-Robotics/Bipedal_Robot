"""The prerequisite must resist stale evidence and boolean-only release."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from design_gate import CHECKS, inspect_gate, require_calculation_release


class GateReleaseTests(unittest.TestCase):
    def setUp(self):
        workspace = Path(__file__).resolve().parents[4]
        scratch = workspace / "tmp"
        scratch.mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix="aarl_design_gate_", dir=scratch)
        self.root = Path(self.temporary.name)
        self.path = self.root / "gate.json"
        artifact = self.root / "verified.txt"
        artifact.write_text("design/pathway fixture", encoding="utf-8")
        digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
        self.gate = {
            "schema": "aarl_all_controller_fsa_design_gate_v1",
            "user_training_approval": {"approved": True, "instruction": "Explicit approval test fixture only"},
            "all_families_required_before_training": True,
            "required_families": ["test_family"],
            "family_evidence": {"test_family": {
                "reviewed_complete": True,
                "checks": {key: {"passed": True, "evidence": "test fixture"} for key in CHECKS},
                "files": [{"path": "verified.txt", "sha256": digest, "role": role}
                          for role in ("design", "pathway_results", "applied_profile", "source")],
            }},
            **{key: True for key in ("training_allowed", "curriculum_allowed",
               "all_family_design_verified", "all_family_scaling_verified",
               "all_family_pathways_verified", "all_family_control_coverage_verified")},
        }

    def tearDown(self):
        self.assertTrue(self.root.resolve().is_relative_to(Path(__file__).resolve().parents[4] / "tmp"))
        self.temporary.cleanup()

    def save(self):
        self.path.write_text(json.dumps(self.gate), encoding="utf-8")

    def test_complete_fixture_releases(self):
        self.save()
        self.assertTrue(require_calculation_release("Hunt_training", self.path)["released"])

    def test_boolean_only_release_is_rejected(self):
        self.gate["family_evidence"] = {}
        self.save()
        with self.assertRaises(RuntimeError):
            require_calculation_release("Hunt_training", self.path)

    def test_completed_design_without_user_approval_stays_closed(self):
        self.gate["user_training_approval"]["approved"] = False
        self.save()
        with self.assertRaises(RuntimeError):
            require_calculation_release("Hunt_training", self.path)

    def test_optional_physics_needs_separate_approval(self):
        self.save()
        for change in ("compliant_tendons", "joint_damping_changes", "muscle_damping_changes"):
            with self.assertRaises(RuntimeError):
                require_calculation_release("curriculum", self.path, physics_changes=[change])
        self.assertTrue(require_calculation_release("curriculum", self.path)["released"])

    def test_changed_source_invalidates_release(self):
        self.save()
        (self.root / "verified.txt").write_text("changed source", encoding="utf-8")
        self.assertFalse(inspect_gate(self.path)["released"])

    def test_missing_control_coverage_is_rejected(self):
        del self.gate["family_evidence"]["test_family"]["checks"]["control_coverage"]
        self.save()
        with self.assertRaises(RuntimeError):
            require_calculation_release("curriculum", self.path)

    def test_one_family_cannot_release_the_other(self):
        self.gate["required_families"].append("second_family")
        self.save()
        self.assertEqual(list(inspect_gate(self.path)["pending_families"]), ["second_family"])

    def test_unknown_operation_rejected(self):
        self.save()
        with self.assertRaises(ValueError):
            require_calculation_release("ignore_gate", self.path)


if __name__ == "__main__":
    unittest.main()
