"""Combine the earlier review table with the historical advisor archive and this revision's readback."""
import html
import json
import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
NOTES = HERE.parent
HISTORY = NOTES/"advisor_feedback_history"
sys.path.insert(0, str(HISTORY))
import review


def text(value):
    return html.escape(str(value), quote=True)


def excerpt(value, limit=240):
    if not value:
        return ""
    if len(value) <= limit:
        return text(value)
    return "<details><summary>"+text(value[:limit])+"…</summary><div>"+text(value)+"</div></details>"


def row(group, identity, location, original, revised, rationale, state):
    fields = (text(group)+"<br><b>"+text(identity)+"</b>", text(location), excerpt(original),
              excerpt(revised), excerpt(rationale), excerpt(state))
    return '<tr data-group="'+text(group)+'">' + "".join("<td>"+f+"</td>" for f in fields) + "</tr>\n"


def anchor_text(value):
    if isinstance(value, dict):
        return str(value.get("text") or value.get("context") or "")
    return str(value or "")


def history_row(s):
    anchor = anchor_text(s.get("anchor"))
    edit = s.get("edit")
    if edit:
        before, after = edit.get("original"), edit.get("proposed")
        fragments = str(edit.get("displayed_fragments",""))
        if before is None or after is None:
            captured = (" Shown fragments: "+fragments) if fragments else ""
            before = "Full original not recovered."+captured if before is None else before
            after = "Full proposal not recovered."+captured if after is None else after
        before = before if before else "(empty original)"
        after = after if after else "(empty proposal)"
    else:
        before, after = anchor or "Exact selected wording unavailable.", s.get("comment","")
    place = f"{s['date']} · {s['author']} · {s['file']}"
    if isinstance(s.get("anchor"),dict) and s["anchor"].get("line"):
        place += f", captured line {s['anchor']['line']}"
    source_status = review.as_text(s.get("evidence_quality")) or "See event source in HISTORICAL_SIGNALS.json."
    return row("Historical advisor",s["id"]+" · "+str(s["score"]),place,before,after,
               s["rationale"],f"Historical {s['type']}. Issue {s['issue_id']}. {source_status}")


def collapse(rows):
    groups = {}
    for s in rows:
        month = s["date"][:7]
        group = groups.setdefault(month,{"events":0,"comments":0,"tracked_edits":0,"deductions":0})
        group["events"]+=1
        group["comments" if s["type"]=="comment" else "tracked_edits"] +=1
        group["deductions"]+=s["score"]
    return dict(sorted(groups.items()))


def main():
    historical = review.read_json(HISTORY/"HISTORICAL_SIGNALS.json")
    rows = review.signals()
    old = historical["signals"]
    decisions = review.read_json(HERE/"PASSAGE_REVIEW.json")
    result = review.make_review(review.DEFAULT_CHAPTERS,rows,decisions)
    assert not result["saved_judgments_stale_or_absent"]
    (HERE/"CHAPTER_SCREEN.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    recent = review.read_json(NOTES/"advisor_scoring_20261118/CLASSIFIER_RESULTS.json")
    # All three inputs are original paragraphs from the previous-day capture.
    # Only pre-November signal text and tags enter the index for this check.
    index = review.EvidenceIndex([r for r in rows if r["date"] < "2026-11-01"])
    expectations = {10:{"mechanism","evidence_and_validation"},12:{"mechanism","evidence_and_validation","result_first"},
                    14:{"gap_and_motivation","contribution_and_credit","reader_use"}}
    retrospective = []
    for p in recent["paragraph_outputs"]:
        if p["file"] != "10-introduction.tex" or p["start"] not in expectations:
            continue
        near = index.nearest(p["input"],limit=3)
        expected = expectations[p["start"]]
        hits = sorted({tag for n in near for tag in n["tags"]} & expected)
        retrospective.append({"line":p["start"],"query":p["input"],"sha256":review.digest(p["input"]),
                              "closest_pre_november_examples":near,"expected_families_from_later_edits":sorted(expected),
                              "families_retrieved":hits})
    assert len(retrospective)==3
    unique = {(r["date"][:7], r["issue_id"] or r["id"]) for r in rows}
    summary={"advisor_signals_by_month":collapse(rows), "count":len(rows),"observed_score":sum(r["score"] for r in rows),
             "scoring_unit":"Each distinct, attributable event appearing in the audited captures, linked by issue within source. This is a recoverable sample.",
             "distinct_issue_month_pairs":len(unique), "screened_units":len(result["passages"]),
             "screened_by_kind":dict(Counter(p["kind"] for p in result["passages"])),
             "dated_local_judgments":len(decisions["passages"]),
             "retrospective_method":"Exact original Introduction text only; all training examples predate November. Match the top three distinct issues against the editorial families observed later. n=3, selected because the advisor edited them. This checks historical retrieval, not classification accuracy, coverage of unmarked prose, or generalization.",
             "november_intro_retrieval":retrospective}

    old_table = (NOTES/"advisor_scoring_20261118/REVIEW_TABLE.html").read_text(encoding="utf-8")
    # Preserve the earlier author directives, failed cases, and other-author review records verbatim.
    old_count = len(re.findall(r'<tr data-group=',old_table))
    new = []
    current_directives = [
        ("NOV18-U1","Does the updated scoring give insight into how to improve your current edits to chapters 6 and 7?",
         "Apply the recurring reader questions to the actual chapter text, retain the engineering scope, and verify the changes.",
         "Ten passages revised. Three already effective passages recorded with reasons for retaining them."),
        ("NOV18-U2","I feel like you should've already had a scoring system in place that looked at the historical edits because sometimes the errors in the way things are presented keep showing up",
         "Recover the prior scoring systems, preserve their original units, and make cross-date evidence part of every future revision.",
         "Added a standing repository instruction, reproducible historical ledger and retrieval, and tests for changed prose and renewed criticism."),
    ]
    for key, request, action, outcome in current_directives:
        new.append(row("Author directive",key,"Ben, November 18 follow-up",request,action,
                       "Ben requested this editorial workflow.",outcome))
    for p in decisions["passages"]:
        new.append(row("Chapter reread",p["id"],f"{p['file']}, line {p['line']}",p["original"],
                       p["proposed_and_applied"],p["finding"]+" Historical examples "+", ".join(p["historical_examples"])+".",
                       f"{p['disposition']}. {p['reviewer_reason']} Local readiness {p['local_readiness']}/4. Advisor acceptance open."))
    new.extend(history_row(s) for s in old)
    headline = (f"<h1>Dissertation feedback and chapter revisions</h1><p>{len(rows)} attributable advisor events in the recovered records "
                f"give an observed score of <strong>−{abs(summary['observed_score'])}</strong> against the original wording. "
                "The November portion is −42. The ten revised paragraphs and three passages retained after comparison have their own local readback. "
                "Select Chapter reread to see exactly what changed. Historical sources preserve imperfect anchors and incomplete edit strings.</p>"
                '<p><a href="REVISION_REVIEW.md">Chapter changes and proof</a> · <a href="PASSAGE_REVIEW.json">Readback and evidence</a> · '
                '<a href="../advisor_feedback_history/HISTORICAL_SIGNALS.json">Historical event sources</a> · '
                '<a href="../advisor_feedback_history/README.md">Scoring rules</a> · <a href="SCORE_SUMMARY.json">Scores and time-separated retrieval check</a> · '
                '<a href="../advisor_scoring_20261118/LIVE_INTRODUCTION_EDITS.json">November 18 Introduction</a> · '
                '<a href="../revision_20261117/LIVE_REVIEW.json">November 17 Review</a>.</p>')
    out = re.sub(r"<h1>.*?</p>\s*<p>.*?</p>(?=\s*<nav>)",lambda _:headline,old_table,count=1,flags=re.S)
    assert headline in out
    out = out.replace("<option>Advisor</option>", "<option>Chapter reread</option><option>Historical advisor</option><option>Advisor</option>",1)
    out = out.replace("<tbody>","<tbody>\n"+"".join(new),1)
    assert len(re.findall(r'<tr data-group=',out)) == old_count+len(new)
    summary["review_table_rows"] = {"prior":old_count,"historical_advisor":len(old),"current_passages":len(decisions["passages"]),
                                     "new_user_directives":len(current_directives),
                                     "total":old_count+len(new)}
    (HERE/"REVIEW_TABLE.html").write_text(out,encoding="utf-8")
    (HERE/"SCORE_SUMMARY.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k!="november_intro_retrieval"},indent=2))
    for test in retrospective:
        print("Intro",test["line"],test["families_retrieved"])


if __name__=="__main__":
    main()
