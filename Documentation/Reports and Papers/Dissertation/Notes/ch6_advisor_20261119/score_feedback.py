"""Extend the dated advisor ledger and check the current passages against the reusable history."""
from collections import Counter
from html import escape
import json
from pathlib import Path
import re
import sys
import tempfile
from zipfile import ZipFile

here = Path(__file__).resolve().parent
notes = here.parent
diss = notes.parent
chapters = diss / "Final/chapters"
sys.path.insert(0, str(notes / "advisor_feedback_history"))
import review

def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

live = review.read_json(here / "LIVE_REVIEW.json")
prior = review.read_json(notes / "advisor_scoring_20261118/CLASSIFIER_RESULTS.json")
dated = notes / "advisor_scoring_20261119"
dated.mkdir(exist_ok=True)

tags = {
    "identify_course_project": ["contribution_and_credit", "reader_context", "ownership_and_narration"],
    "narration": ["ownership_and_narration", "concision_and_repetition"],
    "construction_then_validation": ["narrative_structure", "evidence_and_validation"],
    "connect_problem_to_method": ["gap_and_motivation", "narrative_structure", "concision_and_repetition"],
    "intro_flow": ["gap_and_motivation", "narrative_structure", "mechanism"],
    "tool_as_subject": ["ownership_and_narration", "reader_use"],
}
reasons = {
    "identify_course_project": "Introduce the previous course project once, credit Bolen for the model, and identify Morrow’s collaboration on coupling and tuning.",
    "narration": "Remove the announcement of writing. Retain first person when it identifies who contributed the work.",
    "construction_then_validation": "Give the graphical implementation and the subsequent test of neural–mechanical signal exchange distinct sentences.",
    "connect_problem_to_method": "Connect manual trial cost and the fall on release to the experiment needed to diagnose the cause.",
    "intro_flow": "Connect human and robot mechanics to feedback, the need for controlled comparisons, and the practical limit of the earlier model.",
    "tool_as_subject": "State what the system can simulate. The present-tense verb in the subordinate clause describes an effect within the simulation.",
}
tracked_reasons = [
    "Credit Bolen with building his first model and explicitly credit the work he and Morrow did together.",
    "Use the singular possessive for Bolen’s own first model.",
    "Identify both collaborators for the shared course-project work.",
    "The manual trial cost already establishes the limitation; remove the unsupported comparison of solvers here.",
    "Link the slow runs and hand tuning directly to the resulting difficulty of repeating comparisons.",
    "Make the system the subject of the capability statement.",
    "Use present tense for the effect of a change within the simulation.",
    "Complete the subject–verb agreement in “affects” after the preceding deletion.",
]
signals, others = [], []
for event in live["comments"]:
    row = {k: event[k] for k in ("file","author","time","offset","anchor","text","issue")}
    row.update(id="NOV19-C-"+event["id"], comment_id=event["id"], kind="comment",
               category="advisor_comment", original_score=event["score"], reason=reasons[row["issue"]])
    signals.append(row)
for event in live["edits"]:
    row = event.copy()
    row["kind"] = "replace" if row["original"] and row["proposed"] else "delete" if row["original"] else "insert"
    if row["author"] != "ajh26":
        others.append(row)
        continue
    index = sum(s["category"]=="advisor_tracked_edit" for s in signals)
    row.update(id=f"NOV19-T{index+1}", category="advisor_tracked_edit", original_score=-1,
               anchor=row["proposed"] or row["original"], reason=tracked_reasons[index])
    signals.append(row)
for row in signals:
    row.update(date="2026-11-19", tags=tags[row["issue"]])
for row in others:
    row["date"] = "2026-11-19"
    row["advisor_points"] = 0

all_nov = prior["advisor_signals"] + signals
points = Counter()
for row in all_nov:
    points[row["file"]] += row["original_score"]
assert len(signals) == 15 and len(others) == 3
assert len(all_nov) == len({s["id"] for s in all_nov}) == 56
assert sum(points.values()) == -62
date_score = sum(s["original_score"] for s in signals)
assert date_score == -20
save(dated/"CLASSIFIER_RESULTS.json", dict(
    as_of_date="2026-11-19", method="A human-readable contextual judgment and the standing fixed deductions are applied to each attributable event. Explicit phrase cues and historical word similarity locate passages for rereading. A similarity or cue receives zero advisor points.",
    model_identity="OpenAI Codex, GPT-6 family", rubric=prior["rubric"],
    prompt="Read each advisor message with its original anchor and each tracked replacement in context. Identify the reader’s question, credit the actual author, distinguish writing narration from contributions, and separate current capabilities from proposed experiments. Use the earlier fixed weights on the original prose. Record author edits without advisor points.",
    project=live["url"], source="../ch6_advisor_20261119/LIVE_REVIEW.json",
    historical_ledger="../advisor_scoring_20261118/CLASSIFIER_RESULTS.json",
    advisor_signals=signals, other_authors=others, observed_score_this_capture=date_score,
    observed_events_this_capture=len(signals), november_cumulative_events=len(all_nov),
    observed_total=sum(points.values()), observed_score_by_file=points,
    observed_score_by_date={**prior["observed_score_by_date"],"2026-11-19":date_score}))
with ZipFile(here / "LIVE_SOURCES.zip") as z:
    originals = {name:z.read(name).decode("utf-8") for name in z.namelist()}

# Each judgment identifies actual revised text. The history verifier can recognize it
# after a line move and reopens it if a later advisor comment is anchored there.
cases = [
    ("Course project and ownership", "60-sns.tex", "I built my first neuromechanical model", "identify_course_project", "Model construction is Bolen’s. The two students coupled and tuned the system for a course. The poster has its own attribution."),
    ("Tools and validation in the overview", "60-sns.tex", "I developed the tools in this chapter", "construction_then_validation", "The reader encounters a practical limit, the models and editor, a sensory-to-motion check, and automated experiments in order. Action verbs identify the author’s contributions."),
    ("Why common recordings matter", "60-sns.tex", "The fall could arise", "connect_problem_to_method", "A fall alone leaves three explanations open. Simultaneous recordings can show where the two trials first differ."),
    ("Force at matched poses", "60-sns.tex", "To check whether the converted muscles", "force", "OpenSim and MuJoCo receive the same pose and activation. The comparison reports a large single-pose difference for rectus femoris and a large mean bias for the short head."),
    ("Reason to edit a connectome graphically", "60-sns.tex", "AnimatLab provides a drawing surface", "gui", "The paragraph distinguishes drawing, rendering Python-defined networks, and editing the connection being inspected. The neighboring paragraph gives the scale and an extensor example."),
    ("Library or example builders", "60-sns.tex", "In the graphical editor", "gui", "The editor generates a Simulink circuit. The knee and arm builders supply their separate actuator, sensory, and mechanical blocks."),
    ("Knee coupling test", "60-sns.tex", "In the simulated knee", "reflex", "The electrical stimulus leads through extensor recruitment and force to a measured change in knee angle."),
    ("Arm coupling test", "60-sns.tex", "As the simulated cup filled", "reflex", "A changing external load closes the loop in the other direction. The feedback-on and feedback-off trial use the same loading ramp."),
    ("Journal source for training", "60-sns.tex", "I adapted the four-stage", "hunt", "Bolen adapted the sequence. Hunt, Szczecinski, and Quinn published its four stages and the doubled target period."),
    ("Suspended cycling", "60-sns.tex", "With the same fitted weights", "density", "A measured difference remains, alongside the simultaneous change in drive and feedback and the absence of foot loading."),
    ("Neural timing", "60-sns.tex", "A knee that flexes too little", "density", "A constantly active channel identifies a timing defect in the controller. Separate neural and kinematic recordings identify residual bursts and limited movement in the later suspended trial."),
    ("Supported contact", "60-sns.tex", "On supported level ground", "density", "Under the same assistance and drive, a combined pathway removal disrupted footfall phase while the body continued forward."),
    ("Existing tools", "60-sns.tex", "I developed the SNS-Simscape block library and graphical editor", "conclusion", "The block library, editor, adapted interface, and model preparation are identified with their present operations."),
    ("Existing controlled comparisons", "60-sns.tex", "The calibration and recording scripts", "conclusion", "A researcher can change a circuit or route, hold the other part fixed, and observe coordination, with a scoped timing result."),
    ("Next sensory experiment", "70-overall.tex", "Future researchers can test whether particular", "future", "Apply the existing layer-specific protocol to the six-synergy controller under load; distinguish replay timing from the variability of the moving body."),
    ("Next mechanical experiment", "70-overall.tex", "The force comparison in Appendix", "future", "Vary tendon mechanics or activation separately, then test a characterized BPA force law with recorded pressure."),
    ("Future embedded comparison", "70-overall.tex", "shows how future researchers can test", "future", "A real-time reference and replayed perturbation can test the hardware interface before valve actuation."),
    ("Advisor’s Introduction ending", "10-introduction.tex", "The sensory pathways documented in the database", "intro_flow", "The exact live paragraph connects mechanics and feedback to controlled comparisons, states what the system can simulate, and preserves Ben’s preposition."),
    ("Two muscle discrepancies", "95-AppendixD.tex", "Eight muscles had force-agreement", "force", "The negative scores have an explicit constant-force reference. Separate variation and mean magnitude explain the two cases."),
    ("Controller timing and mechanical outcome", "60-sns.tex", "I changed how the channels", "density", "The timing target was met with prescribed inputs. In a completed suspended run, the third channel produced extra bursts and the rhythm failed its acceptance criteria. The seven-cycle mean knee excursions are observed in the same trial."),
    ("Defense workflow", "95-AppendixD.tex", "follows a saved experiment from model preparation", "gui", "Six stages trace a source model through execution and inspection. The dashed launcher is proposed. The implemented editor draws and builds neural circuits."),
]
extra_tags = {"force":["evidence_and_validation","mechanism"], "gui":["reader_use","gap_and_motivation"],
              "reflex":["evidence_and_validation"], "hunt":["contribution_and_credit"], "density":["result_first","scope_and_uncertainty"],
              "conclusion":["reader_use","contribution_and_credit"], "future":["scope_and_uncertainty","reader_use"]}
rows = review.signals()
index = review.EvidenceIndex(rows)
judgments = []
for i, (title, name, excerpt, issue, rationale) in enumerate(cases, 1):
    hits = [p for p in review.passages(chapters/name) if excerpt in p[3]]
    assert len(hits) == 1, (excerpt, len(hits))
    start, _, _, passage = hits[0]
    similar = index.nearest(passage)[:3]
    judgments.append(dict(id=f"NOV19-READ-{i:02}", title=title, file=name, line=start, issue=issue,
                          tags=tags.get(issue,extra_tags.get(issue,[])), excerpt=excerpt, after=passage,
                          after_sha256=review.digest(passage), historical_examples=similar, judgment=rationale,
                          local_readiness=3, scale="0–4; 4 reserved for author/advisor confirmation"))
save(here/"PASSAGE_REVIEW.json", dict(advisor_event_ids_at_review=[r["id"] for r in rows], passages=judgments))
def diagnostics(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    spans = []
    for label in ("sec:sns_biw2l_results","sec:layer_stimulus_main","sec:sns_ground_results","sec:sns_execution"):
        spans.append(next(i for i, line in enumerate(lines,1) if "\\label{"+label+"}" in line))
    units = [p for p in review.passages(path) if p[2] == "prose"]
    def word_count(text):
        clean = re.sub(r"\\(?:cite\w*|ref|label|path|nolinkurl)\{[^}]*\}", " ", text)
        return len(re.findall(r"\b[A-Za-z]+(?:[’-][A-Za-z]+)*\b",re.sub(r"\\[a-zA-Z]+", " ", clean)))
    return dict(first_person_I=sum(len(re.findall(r"\bI\b",p[3])) for p in units),
        words_by_section={f"6.5.{i+2}": sum(word_count(p[3]) for p in units if lo < p[0] < hi)
                          for i,(lo,hi) in enumerate(zip(spans,spans[1:]))})
with tempfile.TemporaryDirectory() as temporary:
    before = Path(temporary)/"60-sns.tex"
    before.write_text(originals[before.name],encoding="utf-8")
    comparison = dict(live=diagnostics(before), revised=diagnostics(chapters/before.name))
save(here/"PROSE_COUNTS.json", comparison)

def et(s):
    return escape(str(s))
def table(headers, body):
    return '<table><thead><tr>'+"".join("<th>"+et(x)+"</th>" for x in headers)+'</tr></thead><tbody>'+"".join(
        '<tr>'+"".join("<td>"+et(x)+"</td>" for x in row)+'</tr>' for row in body)+'</tbody></table>'
section_names = {"identify_course_project":"Ch. 6 opening and §6.1", "narration":"Ch. 6 opening; §§6.1 and 6.4",
                 "construction_then_validation":"Ch. 6 opening; §6.4.2", "connect_problem_to_method":"§6.1",
                 "intro_flow":"Ch. 1 Motivation, final paragraph", "tool_as_subject":"Ch. 1 Motivation, final sentence"}
def change(row):
    return ('“'+row.get("original","")+'” → “'+row.get("proposed","")+'”')
events_table = table(["File and time (Pacific)","Anchor or tracked change","Advisor comment or rationale","Points","Revised location"],(
    [r["file"]+" · "+r["time"],r["anchor"] if r["kind"]=="comment" else change(r),r.get("text") or r["reason"],
     r["original_score"],section_names[r["issue"]]] for r in signals))
authors_table = table(["File and time","Ben’s change","Advisor points"],([r["file"]+" · "+r["time"],change(r),0] for r in others))
changes_table = table(["Review item","Revised file and line","What the passage now establishes"],(
    [j["title"],j["file"]+" · "+str(j["line"]),j["judgment"]] for j in judgments))
density_table = table(["Section","Earlier narrative words","Revised narrative words"],(
    [label,count,comparison["revised"]["words_by_section"][label]] for label,count in comparison["live"]["words_by_section"].items()))
# The exact anchored live wording and the passage used to address it remain accessible in one table.
case_map = {"identify_course_project":0,"narration":1,"construction_then_validation":1,
            "connect_problem_to_method":2,"intro_flow":17,"tool_as_subject":17}
advice_rows = []
for event in signals:
    j = judgments[case_map[event["issue"]]]
    original = event.get("anchor","") if event["kind"]=="comment" else change(event)
    advice_rows.append([event["id"]+" · "+event["file"]+" · "+event["time"],original,
                       event.get("text") or event["reason"],event["original_score"],
                       j["file"]+" · "+str(j["line"])+" · "+j["judgment"],"Applied locally. Live thread remains with Ben."])
events_table = table(["Advisor event and source","Exact anchor or tracked edit","Comment or purpose","Points","Revised passage and rationale","Disposition"],advice_rows)
other_decisions = [
    ("Ben · first person", "60-sns.tex", "show how I", "Credit concrete contributions, identify collaborative work, and let the tool or experiment state its behavior.", "Chapter 6, opening, §6.1, and conclusion"),
    ("Ben · preserve Introduction", "10-introduction.tex", "the new system can simulate", "The local copy reproduces “combine … with” and the advisor’s ending verbatim.", "Ch. 1, Motivation, final paragraph"),
    ("Ben · preserve Chapter 4 copyedits", "30-torque.tex", "based on bone scans", "Match Ben’s live wording on lines 26 and 28. The subject on line 28 is the singular intersection point. Line 451 matches the full editor capture.", "Ch. 4, lines 26, 28, and 451"),
    ("Introduction comparison · Chapter 6 description", "10-introduction.tex", "graphical robot-model construction", "The chapter demonstrates both operations. CAD supplies the robot geometry, and SNS-Simscape builds the neural circuit. Suggested clause: “and shows how graphical controllers built in SNS-Simscape exchange commands and measurements with a robot model built from CAD.”", "Ch. 1, Dissertation Organization, line 54"),
    ("Ben · neuron terminology", "93-AppendixB.tex", "Graded", "Define nonspiking and spiking before the graphical example. Describe continuous synapses precisely in the hybrid circuit. Update three figure labels and four accessible descriptions.", "§6.4, Appendix B, and Appendix D"),
    ("Ben · placement of two reflex examples", "60-sns.tex", "run_transient.m", "Keep the knee and arm figures. The knee starts with a neural stimulus. The arm starts with a mechanical load. Keep wiring, replay instructions, and the spiking comparison in Appendix D.", "§6.4.2 and Appendix D"),
    ("Ben · graphical programming", "60-sns.tex", "888 neurons", "Use an extensor pathway to illustrate why selecting and changing wiring where it is drawn helps. The SNS-Toolbox renderer displays the Python-defined network. AnimatLab supplies the point of comparison.", "§6.4"),
    ("Ben · shorter results", "60-sns.tex", "iliacus", "Keep the question, matched conditions, result, and material limits in each example. Detailed intervals, scores, and repetitions are in Appendix D.", "§§6.5.2–6.5.4"),
    ("Ben · remove §6.6.3", "60-sns.tex", "I would", "Remove the proposed-experiments subsection. Carry its neural, mechanical, and embedded comparisons into the Chapter 7 future-work account.", "Ch. 7 §7.2.4, Future Work"),
    ("Ben · Chapter 6 conclusion", "60-sns.tex", "a future researcher", "Identify the new block library, graphical editor, preparation and experiment scripts, adapted interface, and present use of controlled comparisons.", "§6.7"),
    ("Ben · Hunt citation", "60-sns.tex", "hunt_inrc_2019", "Use Hunt, Szczecinski, and Quinn (2017), Methods §§2.4.1–2.4.4. Bolen adapted their method to this controller.", "§6.5.1 and Appendix D"),
    ("Ben · physical meaning of force scores", "60-sns.tex", "negative", "Give the force comparison in newtons, show that one muscle disagrees at specific poses while the other has a large mean bias, and state what the constant baseline measures.", "§6.2.1 and Appendix D, Table D2"),
    ("Independent review · sensory scope", "60-sns.tex", "length feedback", "The demonstrated Walker Two-Layer intervention uses contact and force. Qualify the route–stretch claim for circuits that actually receive length feedback.", "§§6.1 and 6.4.3"),
    ("Independent review · measured endpoint", "60-sns.tex", "iliacus", "The phase-specific pulse produced a motor-neuron voltage response. The separate knee experiment establishes force and joint motion.", "§6.5.2"),
    ("Independent review · future protocol", "70-overall.tex", "Future Work", "The reduced-circuit stimulation test already exists. Identify the six-synergy controller under load as the proposed new application.", "Ch. 7, SNS future work"),
    ("Independent review · saved force data", "95-AppendixD.tex", "original framewise", "Locate both 121 × 12 force matrices, match the OpenSim values to the original output, and recompute all twelve statistics. Remove the erroneous archive claim.", "Appendix D, matched-pose test"),
    ("Independent review · definitions", "60-sns.tex", "shared RG", "Identify the rhythm generator’s timing role and gloss Hill-type actuation in the first table.", "§6.2 and Table 6.1"),
    ("Independent review · credit and duplication", "60-sns.tex", "Investigators Meeting", "Identify the poster and its authors. Consolidate the software selector and Hunt attribution in Appendix D.", "§6.1 and Appendix D, air trials"),
    ("Ben · read the source chat and compare both computers", "60-sns.tex", "One syn6 channel", "Read Complete MuJoCo Neural Tasks, verify the exact linked figure, and select the completed suspended latch trial. The directory containing the displayed PNG is named for a different, later experiment. The metrics and raw capture identify the trial actually plotted.", "§6.5.3 and Appendix D"),
    ("Ben · identify what failed", "60-sns.tex", "PF", "The tools exposed a timing error in the original circuit and repeated bursts in a later moving-body test. Report the measured knee motion and the failed rhythm assessment separately. Use the 1.25 m/s Ong simulation and recorded Subject01 walking as different reference sources.", "§6.5.3 and its overlay"),
    ("Ben · defense workflow SVG", "95-AppendixD.tex", "gives the executable sources", "Reflow the six stages for legible print labels. Distinguish the proposed cross-stage launcher from the built circuit editor and model-specific scripts. Keep the existing two-route diagram in Chapter 6.", "Appendix D, execution and reproduction tools"),
]
# Resolve the original lines directly from the prior text. Some concerns were first
# identified in the first local rewrite. These have their exact local-draft wording below.
draft_wording = {
    "Ben · preserve Chapter 4 copyedits":"Earlier local copy: “bones scans” and “The intersection point of the internal crossed links define”. Live: “bone scans” and “The intersection point of the internal crossed links defines”.",
    "Independent review · sensory scope":"Moving an attachment would also change the stretch signal returned to the circuit.",
    "Independent review · measured endpoint":"a command that reaches a muscle",
    "Ben · read the source chat and compare both computers":"One syn6 channel drove its targets throughout the cycle, although it was intended to switch off for part of the step.",
    "Ben · identify what failed":"A circuit triggered by rhythm-generator events passed the timing checks under prescribed sensory replay. When connected to the moving body, it still encountered a knee limit or missed the timing criteria.",
    "Independent review · future protocol":"A second test can stimulate a rhythm-generating, pattern-forming, or motor-neuron pathway separately.",
    "Independent review · saved force data":"The original framewise MuJoCo forces were not archived.",
    "Independent review · credit and duplication":"Hunt supplied the calibration method. […] The software selector for this controller is hunt. The training stages adapt the method of Hunt et al.",
}
decision_rows = []
for title,filename,anchor,after,loc in other_decisions:
    choices=[s.strip() for s in originals[filename].splitlines() if anchor.casefold() in s.casefold()]
    if not choices and title.startswith("Ben · Hunt"):
        choices=[s.strip() for s in originals[filename].splitlines() if "2019" in s and "hunt" in s.lower()]
    earlier=draft_wording.get(title) or (choices[0] if choices else anchor)
    disposition = "Applied locally."
    if "placement" in title:
        disposition = "Both figures retained in Chapter 6 for Ben’s review."
    if title.startswith("Introduction comparison"):
        disposition = "For Ben’s decision. The working Introduction matches his live version."
    decision_rows.append([title,filename+" · "+loc,earlier,after,disposition])
directives_table=table(["Source","Location","Earlier wording or marked passage","Change and rationale","Disposition"],decision_rows)
history_page=(notes/"ch6_first_use_20261118/REVIEW_TABLE.html").read_text(encoding="utf-8")
history_rows={}
for match in re.finditer(r'<tr data-group="([^"]+)">.*?</tr>',history_page,re.S):
    history_rows.setdefault(match[1],[]).append(match[0])
assert len(history_rows["Advisor"])==len(prior["advisor_signals"])==41
previous_columns="<thead><tr>"+"".join("<th>"+c+"</th>" for c in
    ("Source and ID", "Location in the reviewed version", "Earlier wording or anchor", "Proposed wording or action", "Reason and evidence", "Disposition on that date"))+"</tr></thead>"
history_sections="".join("<details><summary>"+escape(group)+" ("+str(len(entries))+")</summary><table>"+
    previous_columns+"<tbody>"+"".join(entries)+"</tbody></table></details>" for group,entries in history_rows.items())
def quantities(text):
    text=re.sub(r"\\(?:cite\w*|path|ref|label|nolinkurl|includegraphics)(?:\[[^\]]+\])?\{[^}]*\}"," ",text)
    text=re.sub(r"\\interlinepenalty=\d+","",text)
    return re.findall(r"(?<![\w.])[-−+]?\d+(?:[.,]\d+)*(?![\w.])",text)
import difflib
numeric_rows = []
for name in ("60-sns.tex","70-overall.tex","93-AppendixB.tex","95-AppendixD.tex"):
    left = originals[name].splitlines()
    right = (chapters/name).read_text(encoding="utf-8").splitlines()
    for action,i,j,k,l in difflib.SequenceMatcher(None,left,right,autojunk=False).get_opcodes():
        if action=="equal":
            continue
        old,new="\n".join(left[i:j]),"\n".join(right[k:l])
        va,vb=quantities(old),quantities(new)
        if va==vb:
            continue
        dest="Same configuration as the live source. See the protocol and interval in the source or the cited Appendix D section."
        if "femoris" in (old+new).lower() or ("mean ratio" in new.lower()):
            dest="Paired OpenSim and MuJoCo forces at the same 121 subject01 poses, 0.5–2.5 s. Source: force_arrays.npz and the original OpenSim static-optimization output. The rectus example is at 1.550 s."
        elif "knee_transient" in (old+new) or "knee_transient" in new.lower():
            dest="Paired knee stimulus and control from the same settled state. See Appendix D, Transient Knee Sensory-Current Response."
        elif "beverage" in (old+new):
            dest="Matched arm load ramps, feedback on and off. Endpoint at 10 s. See Appendix D, Simulink Demonstrations."
        elif "BiW2Lvar" in (old+new) or "app_biw2l_air" in (old+new):
            dest="BiW2Lvar+ on suspended Gait2392. Appendix D retains the paired air trials, loading conditions, and period measurement."
        elif "app_pf_timing" in (old+new):
            dest="Original syn6 nine-cycle PF audit, separate prescribed replay, and completed suspended latch trial. Seven movement cycles from 8.004–13.652 s produce the knee mean-waveform ranges. The eleven ipsilateral rhythm intervals yield the channel counts. See SYN6_OVERLAY.json."
        elif "sec:app_w2l" in (old+new) or "cad_supported" in (old+new):
            dest="Model-specific ground and CAD trials. Appendix D retains the drive, assistance, ledgers, and evaluation windows."
        numeric_rows.append([f"{name} · live {i+1}–{j}, revised {k+1}–{l}",", ".join(va) or "—",", ".join(vb) or "—",dest])
numeric_table=table(["Source block","Numbers before","Numbers after","Configuration, window, or retained location"],numeric_rows)
html = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Chapter 6 · Advisor review and revisions</title><style>body{font:16px/1.55 Georgia,serif;color:#202833;background:#fff;max-width:1150px;margin:3rem auto;padding:0 1.2rem}h1,h2,th{font-family:system-ui,sans-serif}h1{line-height:1.15}h2{font-size:1.25rem;margin-top:2.5rem}table{width:100%;border-collapse:collapse;text-align:left}th,td{padding:.55rem .7rem;border-bottom:1px solid #cbd1d5;vertical-align:top}th{background:#edf2f6}td{overflow-wrap:anywhere}a{color:#15537d}summary{cursor:pointer}@media print{body{margin:0}tr{break-inside:avoid}}</style>
<h1>Chapter 6 advisor review and revisions</h1><p>November 19, 2026 · Dr. Hunt’s Overleaf comments and tracked changes, with Ben Bolen’s edits and directions.</p>
<p>Chapter 6 received six comments and five tracked edits (−15). The Introduction received one comment and three tracked edits (−5). The November 17–19 record contains 56 distinct advisor events totaling −62. The historical prose checker also draws on September and October. Local draft reviews use a separate 0–4 scale. This pass records local verification at 3.</p>
<h2>Advisor comments and tracked changes</h2>"""+events_table+"""
<h2>Author edits</h2>"""+authors_table+"""
<h2>Author directions and independent findings</h2>"""+directives_table+"""
<h2>Revisions and their evidence</h2>"""+changes_table+"""
<h2>Experimental detail in Chapter 6</h2>"""+density_table+"""
<p>Counts cover narrative prose in the three sections. Captions, tabular data, and references are excluded. Appendix D gives the full protocols and measurements.</p>
<details><summary>Numerical changes by source passage</summary>"""+numeric_table+"""</details>
<h2>Earlier comments, decisions, and technical findings</h2>
<p>These records carry forward the original feedback, proposed wording, evidence, and disposition through November 18. Line numbers refer to the draft reviewed on each date. The November 19 rows above give the current decisions for revisited passages.</p>"""+history_sections+"""
<p><a href="LIVE_REVIEW.json">Exact review capture</a> · <a href="../advisor_scoring_20261119/CLASSIFIER_RESULTS.json">Dated scoring ledger</a> · <a href="../advisor_feedback_history/README.md">Standing historical scoring rules</a> · <a href="MANUSCRIPT.diff">Complete line comparison with the live text</a> · <a href="FORCE_CHECK.json">Recalculated muscle forces</a> · <a href="REVISION_REPORT.md">Verification and provenance</a></p></html>"""
(here/"REVIEW_TABLE.html").write_text(html,encoding="utf-8")
print("New signals:",len(signals),date_score,"November:",len(all_nov),sum(points.values()),"Entire history:",len(rows),sum(r["score"] for r in rows))
print(comparison)
