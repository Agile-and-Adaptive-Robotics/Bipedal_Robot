"""Regression checks for historical coverage, reused judgments, and new feedback on old prose."""
import tempfile
import unittest
from pathlib import Path

import review


class HistoryReviewTests(unittest.TestCase):
    def test_historical_coverage_and_newest_dated_deductions(self):
        rows = review.signals()
        self.assertEqual(len(rows), len({r["id"] for r in rows}))
        self.assertTrue(any(r["date"] < "2026-10-01" for r in rows))
        self.assertTrue(any(r["date"].startswith("2026-10") for r in rows))
        captured = review.read_json(review.NOTES/"advisor_scoring_20261118/CLASSIFIER_RESULTS.json")
        ids = {e["id"] for e in captured["advisor_signals"]}
        November = [r for r in rows if r["id"] in ids]
        self.assertEqual((len(November),sum(r["score"] for r in November)),(41,-42))
        newest = [r for r in rows if r["id"].startswith("NOV19-")]
        self.assertEqual((len(newest),sum(r["score"] for r in newest)),(15,-20))
        self.assertTrue(all(r["author"] == "ajh26" for r in newest))

    def test_first_person_is_read_in_context(self):
        words = review.ownership_cues("I show how I built the editor. The new system can simulate a reflex.")
        self.assertIn("roadmap_or_intention", words)
        self.assertIn("I built", words["attributed_contribution"])
        self.assertIn("The new system can", words["named_capability"])
        self.assertNotIn("roadmap_or_intention", review.ownership_cues("I built the editor with Morrow."))

    def test_moved_changed_and_recurring_passage(self):
        body = "The recorded joint torque increased across the measured range of knee flexion."
        path_name = "sample.tex"
        seen = {"file":path_name,"id":"fixture-1","date":"2026-10-01","issue_id":"comparison","score":-2,
                "tags":["evidence_and_validation"],"anchor":"An older torque comparison."}
        saved = {"advisor_event_ids_at_review":[seen["id"]],"passages":[
            {"id":"sample-readback","file":path_name,"after_sha256":review.digest(body)}]}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/path_name
            path.write_text("\n\n"+body+"\n",encoding="utf-8")
            current = review.make_review([path],[seen],saved)
            self.assertEqual(current["passages"][0]["start"],3)
            self.assertEqual(current["passages"][0]["local_judgment"],"sample-readback")
            follow_up = dict(seen,id="fixture-2",date="2026-11-19",anchor=body)
            repeat = review.make_review([path],[seen,follow_up],saved)
            self.assertEqual(repeat["reopened_judgments"],["sample-readback"])
            self.assertIsNone(repeat["passages"][0]["local_judgment"])
            unanchored = review.make_review([path],[seen,dict(follow_up,anchor="")],saved)
            self.assertEqual(unanchored["additional_same_chapter_feedback_to_check"],["fixture-2"])
            path.write_text("The protocol changed.\n",encoding="utf-8")
            stale = review.make_review([path],[seen],saved)
            self.assertEqual(stale["saved_judgments_stale_or_absent"],["sample-readback"])
            self.assertEqual(stale["passages"][0]["status"],"unlabeled")


if __name__ == "__main__":
    unittest.main()
