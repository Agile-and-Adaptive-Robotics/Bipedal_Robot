"""Prune-pattern review from Ben's row examples (2026-10-01, LOCAL ONLY).

Ben's examples (rows 683-822 in his grid) don't map to the export order (his
grid is sorted differently — export row 684 is Brown 1914, which he would
never prune), so the pattern is taken from his CONTENT descriptions instead:
"general neuron modeling … therefore inappropriate for the SADb", "if it
doesn't mention neuroscience", "probably isn't relevant".

Pattern = out-of-scope topical clusters:
  A. General neural-network / ML theory (no biological substrate)
  B. Aerodynamics / flapping-wing engineering
  C. Flow & flight sensing engineering — split by whether the paper shows any
     neuroscience signal (neur*/afferent/reflex/sensorimotor/gangli/CPG...)
  D. Bio-inspired mechanism design (knee joints, leg compliance)
  E. Borderline: fish swimming biomechanics (muscle mechanics, no afferents)
NO Airtable writes — this is a review document with clickable record links.
"""
import json, os, re
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
BASE, TABLE = "appMQTnobUNRytIp7", "tblnnMrZszhboU4uD"
rec_url = lambda rid: f"https://airtable.com/{BASE}/{TABLE}/{rid}"

records = json.load(open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))
layout = json.load(open(os.path.join(HERE, "export", "sadb_layout.json"), encoding="utf-8"))
labels = {int(k): v for k, v in json.load(
    open(os.path.join(HERE, "export", "cluster_labels.json"), encoding="utf-8")).items()}

NEURO = re.compile(r"neur|afferent|reflex|sensorimotor|gangli|cpg|central pattern|"
                   "motor control|electrophysiolog|spike|synap|muscle spindle|"
                   "propriorecept|propriocept|sensill|eme", re.I)

byid = {r["id"]: r for r in records}
def members(cl):
    return sorted((r for r in records if layout.get(r["id"], {}).get("cl", -1) == cl),
                  key=lambda r: (r["primary"] or "?"))

def line(r):
    doi = f" · [doi](https://doi.org/{r['doi']})" if r["doi"] else ""
    return (f"- [{r['primary'] or '?'} {r['year'] or ''} — {r['title'][:88]}]({rec_url(r['id'])})"
            f"{doi}")

out = ["# Prune review — the pattern from your row examples (2026-10-01)",
       "",
       "Your grid row numbers don't match the export order (my row 684 is Brown 1914 —",
       "nothing you'd prune), so I read the pattern from your descriptions instead:",
       "papers whose topic is out of scope for a SENSORY-AFFERENT database — general",
       "neuron/ML modeling, aerospace, non-neural bioinspired engineering. Grouped",
       "below by the named topic clusters, every entry a clickable link to its Airtable",
       "record. **Nothing has been deleted or flagged in Airtable; this is a review",
       "sheet.** Tell me which groups go and I'll flag/execute with your OK (API calls",
       "are scarce until the grace period allows).", ""]

A = [9, 10, 12]      # reservoir computing, neural-net approximation, deep learning
B = [11]             # flapping-wing aerodynamics
C = 7                # bioinspired flow & flight sensing (split by neuro signal)
D = [14, 15]         # knee mechanisms, leg compliance
E = 8                # fish swimming biomechanics

out.append(f"## A. General neural-network / ML theory — {sum(len(members(c)) for c in A)} papers (your rows-811-822 instinct)")
out.append("Reservoir computing, universal-approximator theory, deep-learning methods. No")
out.append("biological substrate, no afferents. Clearest prune group.")
out.append("")
for c in A:
    out.append(f"**{labels[c]}** ({len(members(c))})")
    out += [line(r) for r in members(c)]
out.append("")

out.append(f"## B. Aerodynamics / flapping-wing engineering — {sum(len(members(c)) for c in B)} papers")
out.append("")
for c in B:
    out.append(f"**{labels[c]}** ({len(members(c))})")
    out += [line(r) for r in members(c)]
out.append("")

c7 = members(C)
c7_neuro = [r for r in c7 if NEURO.search((r["title"] or "") + " " + (r.get("notes") or ""))]
c7_plain = [r for r in c7 if r not in c7_neuro]
out.append(f"## C. Bioinspired flow & flight sensing — {len(c7)} papers, split by neuroscience signal (your 686/708 'if it doesn't mention neuroscience')")
out.append("")
out.append(f"### C1. No neuroscience signal in title+notes — {len(c7_plain)} prune-leaning")
out += [line(r) for r in c7_plain]
out.append("")
out.append(f"### C2. Mentions neural/sensorimotor concepts — {len(c7_neuro)} keep-leaning (your call)")
out += [line(r) for r in c7_neuro]
out.append("")

out.append(f"## D. Bio-inspired mechanism design — {sum(len(members(c)) for c in D)} papers (matches your 'probably isn't relevant')")
out.append("")
for c in D:
    out.append(f"**{labels[c]}** ({len(members(c))})")
    out += [line(r) for r in members(c)]
out.append("")

e8 = members(E)
out.append(f"## E. Borderline: fish swimming biomechanics — {len(e8)} papers")
out.append("Muscle power/work kinematics (bluegill, cod). No afferent focus, but muscle")
out.append("mechanics borders on the SADb's modeling use. Listed separately — NOT")
out.append("assumed pruneable.")
out += [line(r) for r in e8]
out.append("")

nA = sum(len(members(c)) for c in A)
nB = sum(len(members(c)) for c in B)
nD = sum(len(members(c)) for c in D)
tot_strict = nA + nB + len(c7_plain) + nD
out.append(f"## Counts")
out.append(f"- Strict pattern (A + B + C1 + D): **{tot_strict} papers**")
out.append(f"- Plus borderline E: {tot_strict + len(e8)}")
out.append(f"- Combined with the 38 already flagged (isolated + uncurated + uncited),")
out.append(f"  deleting everything on the table would free up to **{tot_strict + len(e8) + 38}** of the")
out.append(f"  ~390 needed to reach the 1000-record cap (base grew to ~1389 with the twins).")
out.append(f"- Reminder: the earlier 38 Prune-Status flags stay as-is in Airtable — no action taken per your word.")

path = os.path.join(HERE, "prune_review_20261001.md")
open(path, "w", encoding="utf-8").write("\n".join(out) + "\n")
print(f"wrote {path}: strict {tot_strict} (+{len(e8)} borderline, +38 prior flags)")
