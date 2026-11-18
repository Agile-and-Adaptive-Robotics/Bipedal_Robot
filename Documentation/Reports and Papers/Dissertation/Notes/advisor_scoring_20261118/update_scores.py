"""Add the observed Introduction edits to the November 17 classifier ledger and review table."""
from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
from html import escape
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import zipfile

here = Path(__file__).resolve().parent
previous = here.parent / "revision_20261117"
intro = "10-introduction.tex"
capture_path = here / "LIVE_INTRODUCTION_EDITS.json"
capture = json.loads(capture_path.read_text(encoding="utf-8"))
prior_path = previous / "CLASSIFIER_RESULTS.json"
prior = json.loads(prior_path.read_text(encoding="utf-8"))
assert prior["observed_total"] == -6 and len(prior["advisor_signals"]) == 5

# The individual edit and the interpreted editorial issue are different units.
# These six judgments are the language model's classifications of the original
# passages after reading both wordings. The score of each edit follows Ben's rubric.
issue_groups = [
    ("intro_actuator_model", 1, 5, "Resting length and the force experiment",
     "Identify the relationship missing from the empirical model, the data collected on the jig, and the force model those data enabled."),
    ("intro_torque_mechanism", 6, 11, "How artificial muscles move a limb",
     "Introduce the joint moment and its role in movement before describing how force reaches the skeleton and how the attachments affect it."),
    ("intro_torque_prediction", 12, 16, "Why rigid-body predictions miss measured knee behavior",
     "Relate compliance and bending to the torque prediction, then describe the shortfall observed on the knee rig."),
    ("intro_route_design", 17, 21, "What the knee tests and redesigned routes achieved",
     "State what the model predicts and identify the meet-or-exceed human torque criterion for the redesigned routes."),
    ("intro_neural_motivation", 22, 29, "Why the controller project needs a literature collection",
     "Connect the neural-controller goal to the need to gather evidence about sensory control across species, preparations, and protocols."),
    ("intro_sadb_contribution", 30, 36, "Ownership and use of the SADb",
     "Name the author as the database developer and describe how readers can organize and explore the literature."),
]
reasons = [
    "Identifies the resting-length relationship missing from the previous force models.",
    "Adds the purpose of collecting isometric measurements.",
    "Identifies the experimental data as the basis for model development.",
    "Specifies which force the model predicts.",
    "Tells the reader what the muscle chapter presents.",
    "Shortens the opening of the assembled-joint problem.",
    "Explains how muscle force produces a joint moment and why the moment matters for movement.",
    "Connects artificial-muscle force to transmission through brackets and tendons.",
    "Identifies the structure that contacts the actuator.",
    "Removes a second consequence from the list of physical mechanisms.",
    "States the consequence for force and torque on the limb.",
    "Refers directly to the compliance mechanisms just described.",
    "Removes the punctuation left by restructuring the list.",
    "Specifies that the rigid-body comparison concerns torque in robot models.",
    "Contrasts the prediction with the torque the assembled joint produces.",
    "Describes an observed knee-test shortfall in place of a hypothetical possibility.",
    "Explains the mechanical purpose of the model developed with the test jig.",
    "Identifies the redesigned routes as artificial musculotendon routes.",
    "States the design criterion as meeting or exceeding the human comparison.",
    "Names the human muscles as the comparison for the redesigned routes.",
    "Introduces the chapter reference as the location of the design procedure.",
    "Introduces the human-like neural controller as a goal of the bipedal project.",
    "Links the evidence collection to the controller the project aims to build.",
    "Uses the present tense for the continuing need for this information.",
    "Describes the requested literature resource in terms a new reader can follow.",
    "Includes nervous-system organization as well as the effects of sensory input.",
    "Removes an early list of source types to keep the need for the resource in focus.",
    "Makes the gap in the collected understanding of neural control explicit.",
    "Shortens the evidence-gap sentence after naming its experimental contexts.",
    "Attributes the development of the SADb directly to the author.",
    "Removes the redundant passive verb after the contribution is attributed.",
    "Reframes the database as a means to organize the literature.",
    "Completes the verb form required by the rephrased purpose.",
    "Names the literature as the material being organized.",
    "Makes the chapter itself the subject of the final sentence.",
    "Describes the construction, organization, and navigation features readers will encounter.",
]

def sha(data):
    return hashlib.sha256(data).hexdigest()

def ui_checksum(text):
    data = text.encode("utf-16-le")
    value = 2166136261
    for i in range(0, len(data), 2):
        value = ((value ^ int.from_bytes(data[i:i+2], "little")) * 16777619) & 0xffffffff
    return f"{value:08x}", len(data) // 2

def json_out(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

events = capture["entries"]
assert len(events) == len(reasons) == 36
assert len({(e["offset"],e["kind"],e["original"],e["proposed"]) for e in events}) == len(events)
assert Counter(e["kind"] for e in events) == {"replace":21, "insert":9, "delete":6}
assert {e["file"] for e in events} == {intro} and {e["author"] for e in events} == {"ajh26"}
assert not any(s["file"] == intro for s in prior["advisor_signals"])
card_text = [[e[k] for k in ("offset","kind","time","original","proposed")] for e in events]
compact = json.dumps(card_text, ensure_ascii=False, separators=(",",":"))
assert ui_checksum(compact)[0] == capture["ui_integrity"]["cards"]
with zipfile.ZipFile(previous / "LIVE_SOURCE_TEXT.zip") as baseline_zip:
    original = baseline_zip.read(intro).decode("utf-8-sig").replace("\r\n","\n")
baseline_lines = original.splitlines()
proposed = original
group_by_id = {}
for key, first, last, _, _ in issue_groups:
    for i in range(first, last+1):
        assert i not in group_by_id
        group_by_id[i] = key
scores = defaultdict(list)
new_signals = []
for i, e in enumerate(events, 1):
    offset = e["offset"]
    old, new = e["original"], e["proposed"]
    assert proposed[offset:offset+len(old)] == old, e["id"]
    line = proposed.count("\n", 0, offset) + 1
    proposed = proposed[:offset] + new + proposed[offset+len(old):]
    scores[line].append(e["id"])
    new_signals.append({**e,"category":"advisor_tracked_edit","original_score":-1,
                        "issue":group_by_id[i],"reason":reasons[i-1],"live_line":line})
live_lines = proposed.splitlines()
assert {n:len(ids) for n,ids in scores.items()} == {10:5,12:16,14:15}
for n, control in capture["ui_integrity"]["proposed_source_lines"].items():
    assert ui_checksum(live_lines[int(n)-1]) == (control["hash"], control["length"])
for s in new_signals:
    s["live_tex_excerpt"] = live_lines[s["live_line"]-1]
    s["related_issue_grouping"] = "Interpretation of related edits by the language-model classifier."

# Keep the earlier predictions, including unmarked paragraphs, intact. Chapter 1
# contributes exactly the three observations covered by this scoring request.
result = deepcopy(prior)
result["model_runs"] = {"2026-11-17":prior["model_identity"],
                        "2026-11-18":"OpenAI Codex, GPT-6 family, main session. Exact checkpoint and sampling settings unavailable."}
result["model_identity"] = "OpenAI Codex, GPT-6 family. The model_runs field identifies the reviewer for each dated batch."
result["method_20261117"] = prior["method"]
result["method"] = ("The November 17 classifier and its outputs are retained. On November 18 the main-session language model read the original "
                    "passages, the 36 expanded changes, and the proposed passages in context. It assigned the issue group and reason for each "
                    "observed edit. This script applies the fixed scoring rule, replays the text against its recorded offsets, and calculates the totals.")
result["prompt_20261117"] = prior["prompt"]
result["prompt"] = ("For each distinct change by ajh26 in the November 18 Introduction review, record a score of -1 on the original draft. "
                    "Count an insertion, deletion, or replacement as one event, including punctuation and word fragments. For each event identify "
                    "its editorial purpose from the prior and proposed sentences. Group related changes by the reader's question they address. "
                    "Carry the November 17 observed scores forward once. Keep model predictions and other authors' edits outside the observed total. "
                    "Record the exact original and proposed wording. Preserve the advisor's own wording and leave publication decisions with Ben.")
result["captures"] = [
    {"date":"2026-11-17","file":"../revision_20261117/LIVE_REVIEW.json","sha256":prior["review_capture_sha256"]},
    {"date":capture["date"],"file":capture_path.name,"sha256":sha(capture_path.read_bytes())}
]
result["capture_utc_20261117"] = result.pop("capture_utc")
result["capture_date_20261118"] = capture["date"]
result["as_of_date"] = capture["date"]
result["review_capture"] = capture_path.name
result["review_capture_sha256"] = sha(capture_path.read_bytes())
result["known_anchor_rule_20261117"] = prior["known_anchor_rule"]
result["known_anchor_rule"] = ("The earlier run covers Chapters 6 and 7 and Appendices A–E. The November 18 Introduction observations are "
                               "on lines 10, 12, and 14. All 36 edits replay against the November 17 text, and the three resulting paragraphs "
                               "match the text in the open editor. Other Introduction passages were outside this additional scoring pass.")
result["source_snapshots"][intro] = {
    "original_capture":"../revision_20261117/LIVE_SOURCE_TEXT.zip::10-introduction.tex",
    "original_utf8_lf_sha256":sha(original.encode()),
    "verified_proposed_paragraphs":{str(n):sha(live_lines[n-1].encode()) for n in scores}}
result["advisor_signals"].extend(new_signals)
result["observed_total"] = sum(s["original_score"] for s in result["advisor_signals"])
assert result["observed_total"] == -42
result["observed_score_by_file"] = {intro:-36,"60-sns.tex":-6}
result["observed_score_by_date"] = {"2026-11-17":-6,"2026-11-18":-36}
result["source_20261117"] = prior["source"]
result["source"] = "November 17 review overview and archived source text, plus the November 18 Introduction review cards and three live editor paragraphs."
for key, first, last, title, purpose in issue_groups:
    result["issue_catalog"][key] = title + ". " + purpose
    result["observed_unique_issues"][key] = [e["id"] for e in events[first-1:last]]
result["prose_counts"]["observed_tracked_edit_anchor"] += len(scores)
result["chapter_unit_counts"][intro] = len(scores)
for n, ids in scores.items():
    result["paragraph_outputs"].append({
        "file":intro,"start":n,"end":n,"section":"Motivation","input":baseline_lines[n-1],
        "input_sha256_utf8_lf":sha(baseline_lines[n-1].encode()),"advisor_label":None,"predicted_severity":None,
        "predicted_issue":None,"model_reason":None,"status":"observed_tracked_edit_anchor",
        "advisor_signals":ids,"observed_score_sum":-len(ids),
        "recorded_live_proposed_paragraph":live_lines[n-1],"git_start":n})
assert result["paragraph_outputs"][:268] == prior["paragraph_outputs"]
assert len(result["advisor_signals"]) == 41
assert sum(result["prose_counts"].values()) == sum(result["chapter_unit_counts"].values()) == len(result["paragraph_outputs"]) == 271
json_out(here/"CLASSIFIER_RESULTS.json",result)

# Extend the existing single-table review, retaining its rows, search, and filters.
html = (previous/"REVIEW_TABLE.html").read_text(encoding="utf-8")
html = html.replace("<h1>Dissertation revision review</h1>", "<h1>Dissertation revision review through November 18</h1>")
html = html.replace("The advisor gave five scored signals on two issues, totaling −6.",
                    "The cumulative advisor score is <strong>−42</strong>. Chapter 1 contributed 36 tracked edits (−36). Chapter 6 contributed three tracked edits and two comment messages (−6).")
html = html.replace("5 advisor signals,", "41 advisor signals,")
html = html.replace('href="LIVE_REVIEW.json"', 'href="../revision_20261117/LIVE_REVIEW.json"')
html = html.replace('>Verbatim live capture</a>', '>November 17 review</a> · <a href="LIVE_INTRODUCTION_EDITS.json">November 18 Introduction edits</a> · <a href="CLASSIFIER_RESULTS.json">Cumulative scores</a>')
for name in ("INDEPENDENT_REVIEW.md","TECHNICAL_READER.md"):
    html = html.replace(f'href="{name}"',f'href="../revision_20261117/{name}"')

def current_disposition(markup, issue_id, value):
    pattern = r'(<b>'+re.escape(issue_id)+r'</b></td>(?:(?!</tr>).)*<td>)[^<]*(</td></tr>)'
    updated, hits = re.subn(pattern,lambda m:m[1]+escape(value)+m[2],markup,flags=re.S)
    assert hits == 1
    return updated

html = current_disposition(html,"IR18","Ben corrected both typos in Overleaf. Verified November 18.")
html = current_disposition(html,"X01","Ben retained line 451. Its opening “At each knee angle” supports the pose-dependent reading in Appendix A.")
def wording(s):
    return escape(s) if s else '<span class="empty">Empty string</span>'
rows=[]
for s in new_signals:
    title = next(g[3] for g in issue_groups if g[0]==s["issue"])
    cells=[f'Advisor<br><b>{s["id"]} · −1</b>',
           f'10-introduction.tex, line {s["live_line"]}, offset {s["offset"]}<br>{escape(s["time"])} Pacific · ajh26',
           wording(s["original"]),wording(s["proposed"]),escape(s["reason"])+" "+escape(title)+".",
           "Scored against the original. Proposed wording is visible in the Overleaf Review panel."]
    rows.append('<tr data-group="Advisor">'+"".join("<td>"+v+"</td>" for v in cells)+"</tr>")
html = html.replace("<tbody>","<tbody>\n"+"\n".join(rows)+"\n",1)

class TableCheck(HTMLParser):
    def __init__(self):
        super().__init__()
        self.groups=Counter()
        self.links=[]
        self.tables=0
    def handle_starttag(self, tag, attrs):
        a=dict(attrs)
        if tag=="table": self.tables+=1
        if tag=="tr" and "data-group" in a: self.groups[a["data-group"]]+=1
        if tag=="a" and "href" in a: self.links.append(a["href"])

old_check=TableCheck()
old_check.feed((previous/"REVIEW_TABLE.html").read_text(encoding="utf-8"))
check=TableCheck()
check.feed(html)
assert check.tables==1
assert check.groups == old_check.groups + Counter({"Advisor":36})
assert all((here/f).exists() for f in check.links)
assert all(html.count(e["id"]+" · −1")==1 for e in events)
(here/"REVIEW_TABLE.html").write_text(html, encoding="utf-8")
print(f"{len(events)} new edits; {len(result['advisor_signals'])} cumulative signals; {result['observed_total']} points; {sum(check.groups.values())} review rows.")
