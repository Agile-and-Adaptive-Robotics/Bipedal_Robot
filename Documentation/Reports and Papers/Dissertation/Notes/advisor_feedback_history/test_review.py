"""Regression checks for historical coverage, reused judgments, and new feedback on old prose."""
import tempfile
import unittest
from pathlib import Path

import review


class HistoryReviewTests(unittest.TestCase):
    def test_advisor_style_review_resumes_after_the_author_reserved_ethics_section(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "15-background.tex"
            path.write_text(
                "\\section{Methods}\n\nMeasure the movement.\n\n"
                "\\section{Ethics}\\label{sec:bg_ethics}\n\nThe model replaces the animal experiment.\n\n"
                "\\section{Next section}\n\nMeasure the force.\n", encoding="utf-8")
            units = review.make_review([path], review.signals())["passages"]
        self.assertEqual([unit["status"] for unit in units],
                         ["unlabeled", "author-reserved ethical argument", "unlabeled"])
        self.assertEqual(units[1]["historical_examples"], [])
        self.assertIsNone(units[1]["local_judgment"])

    def test_results_openings_flag_inventory_and_instructions_but_keep_quantitative_findings(self):
        text = r"""\section{Results}

The completed Explorer contains 1,393 records. Its 943 \emph{Papers} are bibliographic entries.

\subsection{From a Pathway Record to a Controller Rule}

Does input from one ankle extensor reinforce activity in another? A reader can repeat the lookup by opening the file.

\subsection{Torque comparison}

The measured torque exceeded the target at 15 configurations. The measurements give the margin for redesign.

\section{Methods}

Open the file and choose the pathway.
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"chapter.tex"
            path.write_text(text, encoding="utf-8")
            tested = review.make_review([path], review.signals())["results_opening_review"]
        self.assertEqual([p["cues"] for p in tested["openings"]],
                         [["inventory leads the result"], ["operating instructions lead the result"], []])
        self.assertEqual({r["type"] for r in tested["basis"]}, {"comment", "tracked_edit"})
        tracked = next(r for r in tested["basis"] if r["id"] == "INTRO-T36")
        self.assertIn("sort and understand", tracked["edit"]["proposed"])

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
