"""Round-2 prune sweep + duplicate scan (2026-10-01, LOCAL ONLY).

Ben's new patterns (content, not row numbers — those don't map):
  1. CS-style "neural networks" (prediction/classification/optimization) that
     are NOT biologically inspired
  2. Traditional control theory (PID/MPC/robust/adaptive controller) with no
     neural/pathway basis
  3. General neuron models (spiking-neuron math etc.) without biology
Plus: duplicate scan — same-author pairs with near-identical titles or same
author+year (his rows 275/276 hint). Excludes records already deleted.
Output: prune_review_round2.md with clickable links; nothing is executed.
"""
import json, os, re
from collections import defaultdict
from difflib import SequenceMatcher

HERE = os.path.dirname(os.path.abspath(__file__))
BASE, TABLE = "appMQTnobUNRytIp7", "tblnnMrZszhboU4uD"
rec_url = lambda rid: f"https://airtable.com/{BASE}/{TABLE}/{rid}"

records = json.load(open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))
deleted = {r["id"] for r in json.load(open(os.path.join(HERE, "deleted_records_20261001.json"),
                                           encoding="utf-8"))}
records = [r for r in records if r["id"] not in deleted]
print("live papers:", len(records))

BIO = re.compile(r"locomot|walk|gait|cpg|central pattern|reflex|afferent|spinal|neuro|insect|"
                 r"muscle|motor|sens|biolog|bio-?inspired|neuromechan|animal|frog|fish|lamprey|"
                 r"salamander|stick insect|cockroach|crayfish|locust|moth|limb|legged|posture|"
                 r"swim|scratch|vision guided by|optomotor|antennal|mantis|hexapod|rat |mouse|cat ",
                 re.I)
CS_NN = re.compile(r"neural network|deep learning|machine learning|reinforcement|backpropag|"
                   r"lstm|recurrent network|echo state|reservoir|convolutional|radial basis|"
                   r"extreme learning|spiking neuron|artificial neuron|neuron model", re.I)
CTRL = re.compile(r"\bpid\b|model predictive|robust control|adaptive control|optimal control|"
                  r"trajectory tracking|lqr|h[ _]infinity|sliding mode|feedback control|"
                  r"controller design|gain scheduling|impedance control", re.I)

def text(r):
    return (r["title"] or "") + " " + (r.get("notes") or "")[:400]

cs_nn, ctrl, neuron = [], [], []
for r in records:
    t = text(r)
    if CS_NN.search(t) and not BIO.search(t):
        cs_nn.append(r)
    elif CTRL.search(t) and not BIO.search(t) and not CS_NN.search(t):
        ctrl.append(r)
    elif re.search(r"neuron model|neuronal model|spiking model|izhikevich|hodgkin|huxley|"
                   r"fitzhugh|hindmarsh", t, re.I) and not BIO.search(t):
        neuron.append(r)

# duplicate scan: same primary author, title similarity or same year
pairs = []
by_author = defaultdict(list)
for r in records:
    a = (r["primary"] or "").strip().lower()
    if a:
        by_author[a].append(r)
for a, rs in by_author.items():
    for i in range(len(rs)):
        for j in range(i + 1, len(rs)):
            t1, t2 = (rs[i]["title"] or "").lower(), (rs[j]["title"] or "").lower()
            sim = SequenceMatcher(None, t1, t2).ratio()
            same_year = rs[i]["year"] == rs[j]["year"] and rs[i]["year"]
            if sim > 0.72 or (sim > 0.5 and same_year):
                pairs.append((rs[i], rs[j], round(sim, 2), same_year))

def line(r):
    doi = f" · [doi](https://doi.org/{r['doi']})" if r["doi"] else ""
    return f"- [{r['primary'] or '?'} {r['year'] or ''} — {r['title'][:88]}]({rec_url(r['id'])}){doi}"

out = ["# Prune review — ROUND 2 (Ben's 2026-10-01 patterns) — nothing executed",
       "",
       "Patterns from your message, applied to the 855 live papers (row numbers",
       "in the grid don't map to records, so these are content matches on",
       "title + curation note; each entry links to its Airtable record).", ""]

def section(title, items, note):
    out.append(f"## {title} — {len(items)} papers")
    out.extend(["", note, ""])
    out.extend([line(r) for r in sorted(items, key=lambda r: (r["primary"] or "?"))] or ["- (none matched)"])
    out.append("")

section("CS-style neural networks (not biologically inspired)", cs_nn,
        "Prediction/classification/optimization NN usage with no biological framing — "
        "your rows-736-759 / 'computer-science neural networks' pattern.")
section("Traditional control (no neural basis)", ctrl,
        "PID/MPC/robust/adaptive controller papers — your rows 844/854/827/843 pattern. "
        "NOTE: Park 2003 and Wang 2020 are human postural-control EXPERIMENTS (biological "
        "subjects, EMG/force data) — auto-matched on 'control' wording but likely IN scope; "
        "Hyun 2014 is a classical hierarchical controller using proprioceptive signals — "
        "the robotics-boundary case you described. Giesseler and Wang 2013 are pure "
        "aerospace control — clear prunes.")
section("General neuron models (no biology)", neuron,
        "Neuron-modeling math without a biological/preparation context.")

out.append(f"## Duplicate candidates — {len(pairs)} pairs (your rows-275/276 hint)")
out.append("")
out.append("Same primary author + similar title (or same year). AUTO-VERDICT: a bioRxiv "
           "DOI (10.1101/...) on one side = true preprint/published duplicate (prune the "
           "preprint); titles differing by part numbering (I/II), hindlimb/forelimb, or "
           "data/model = companion papers, KEEP BOTH.")
out.append("")
for a, b, sim, sy in sorted(pairs, key=lambda p: -p[2]):
    doia, doib = a["doi"], b["doi"]
    if doia.startswith("10.1101/") or doib.startswith("10.1101/"):
        verdict = "**TRUE DUPLICATE (preprint vs published — prune the 10.1101 side)**"
    elif re.search(r"\b(i{1,3})\b.*\b(i{1,3})\b", a["title"] + "|" + b["title"]) or \
         "hindlimb" in a["title"].lower() and "forelimb" in b["title"].lower() or \
         "forelimb" in a["title"].lower() and "hindlimb" in b["title"].lower():
        verdict = "companion papers (different studies — keep both)"
    else:
        verdict = "review: companion or duplicate?"
    out.append(f"- sim {sim}{f' · same year {sy}' if sy else ''} — {verdict}")
    out.append(f"  - [{a['primary']} {a['year']} — {a['title'][:80]}]({rec_url(a['id'])}) ({doia[:44]})")
    out.append(f"  - [{b['primary']} {b['year']} — {b['title'][:80]}]({rec_url(b['id'])}) ({doib[:44]})")
out.append("")

path = os.path.join(HERE, "prune_review_round2.md")
open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")
print(f"wrote {path}: cs_nn {len(cs_nn)} | ctrl {len(ctrl)} | neuron {len(neuron)} | dup pairs {len(pairs)}")
