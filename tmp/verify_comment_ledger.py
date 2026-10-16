"""Verify comment_ledger.json: JSON parse, structural rules, quote-to-file placement.
Run: C:\\Users\\Ben Bolen\\.conda\\envs\\myo\\python.exe tmp\\verify_comment_ledger.py
"""
import json, os, re, sys

ROOT = r"D:\Github\Bipedal_Robot"
LEDGER = os.path.join(ROOT, r"Documentation\Reports and Papers\Dissertation\Notes\revision_campaign_20261015\comment_ledger.json")
CH = os.path.join(ROOT, r"Documentation\Reports and Papers\Dissertation\ProofFinal\chapters")
MAIN = os.path.join(ROOT, r"Documentation\Reports and Papers\Dissertation\ProofFinal\main.tex")

files = {}
for fn in os.listdir(CH):
    if fn.endswith(".tex"):
        files[fn] = open(os.path.join(CH, fn), encoding="utf-8", errors="replace").read()
files["main.tex"] = open(MAIN, encoding="utf-8", errors="replace").read()

def norm(s):
    return re.sub(r"\s+", " ", s).strip()

# ---------- 1. parse ----------
with open(LEDGER, encoding="utf-8") as f:
    ledger = json.load(f)
print("JSON PARSE OK — %d entries" % len(ledger))

# ---------- 2. structure ----------
expected_ids = ["C%02d" % i for i in range(1, 47)] + ["T%02d" % i for i in range(1, 16)]
ids = [e["id"] for e in ledger]
problems = []
if ids != expected_ids:
    problems.append("id list/order mismatch: %s" % (set(expected_ids) ^ set(ids)))

ACTIONS = {"rewrite", "expand", "move", "cut", "fix-math", "accept-TC", "layout"}
HOMES = {"intro", "background", "muscles", "torque", "ankle", "database", "sns", "overall", "backmatter"}
AUTHORS = {"ajh26", "mcneal", "bbolen83"}
BASE_FIELDS = {"id", "author", "quotedHighlight", "oldSection", "action", "weight", "newHome", "anchorConfidence"}
STRUCT3 = {"C12", "C13", "C14", "C15", "C22", "C23", "C30"}
DETACHED = {"C10", "C11"} | {"C%02d" % i for i in range(24, 38)} | {"C%02d" % i for i in range(38, 44)} | {"C45", "C46"}
GROUPED_TC = {"T01", "T08", "T09", "T10", "T11", "T12", "T13", "T14", "T15"}

for e in ledger:
    i = e["id"]
    extra = set(e) - BASE_FIELDS
    if extra - {"note"}:
        problems.append("%s: unexpected fields %s" % (i, extra))
    for fld in BASE_FIELDS:
        if fld not in e:
            problems.append("%s: missing field %s" % (i, fld))
    if e.get("action") not in ACTIONS: problems.append("%s: bad action %r" % (i, e.get("action")))
    if e.get("newHome") not in HOMES: problems.append("%s: bad newHome %r" % (i, e.get("newHome")))
    if e.get("author") not in AUTHORS: problems.append("%s: bad author %r" % (i, e.get("author")))
    if e.get("anchorConfidence") not in {"quoted", "inferred"}: problems.append("%s: bad anchorConfidence" % i)
    if e.get("weight") not in (1, 2, 3): problems.append("%s: bad weight %r" % (i, e.get("weight")))
    # weight rules
    if i in STRUCT3 and e["weight"] != 3: problems.append("%s: structural comment must be weight 3" % i)
    if e["author"] == "mcneal" and e["weight"] != 3: problems.append("%s: mcneal must be weight 3" % i)
    if e["author"] == "bbolen83" and e["weight"] != 1: problems.append("%s: bbolen83 must be weight 1" % i)
    # anchor rules
    want = "inferred" if (i in DETACHED or i in GROUPED_TC) else "quoted"
    if e["anchorConfidence"] != want:
        problems.append("%s: anchorConfidence %s, expected %s" % (i, e["anchorConfidence"], want))
    if i in DETACHED and e["quotedHighlight"] != "":
        problems.append("%s: detached comment must have empty quotedHighlight" % i)

print("STRUCTURE: %d problem(s)" % len(problems))
for p in problems: print("  " + p)

# ---------- 3. quote placement ----------
print("\nQUOTE PLACEMENT (normalized-whitespace substring match in ProofFinal):")
misses = []
for e in ledger:
    i, q = e["id"], e["quotedHighlight"]
    if not q or i in DETACHED or i in GROUPED_TC or i in {"T04", "T05"}:
        continue  # empty, detached, grouped concatenations, or added-text (position not in source)
    fname = e["oldSection"].split(" \u2192 ")[0]
    nq = norm(q)
    hit = nq in norm(files.get(fname, ""))
    where = fname if hit else "NOWHERE"
    if not hit:
        for fn, txt in files.items():
            if nq in norm(txt):
                where = fn + " (other file)"
                hit = True
                break
    print("  %-4s %-28s %s" % (i, where, "OK" if hit else "MISS"))
    if not hit: misses.append(i)

print("\nGROUPED-TC FRAGMENT SPOT CHECKS (concatenated deletions; fragments may be reworded in current tree):")
frag = {
 "T01": ("10-introduction.tex", ["\\section{Motivation}", "three gaps stand between the existing robot design", "The present work builds directly on a series of studies by the author and collaborators."]),
 "T08": ("20-methods.tex", ["The bracket stiffness matrix is modeled as a diagonal matrix", "bktFrame.pdf", "The goodness-of-fit metrics used to judge the optimization results"]),
 "T09": ("30-results.tex", ["The AnimatLab phase-1 run described in Section", "The identification work outlined at the close of the experimental chapters", "The route redesign built on these values"]),
 "T10": ("40-discussion.tex", ["The reason for BPA models predicting higher force", "biological muscle optimum fiber length", "rigid body mechanics simplifications"]),
 "T11": ("50-futurework.tex", ["regularize the behavior of different", "proprioceptive feedback", "Jetson"]),
 "T12": ("60-conclusion.tex", ["select actuator diameter and resting length, and to evaluate", "is realizable in hardware, and identifying"]),
 "T13": ("92-AppendixA.tex", ["MATLAB code", "compressive"]),
 "T14": ("94-AppendixC.tex", ["Regenerating the results files updates this table."]),
 "T15": ("main.tex", ["Creating"]),
}
fragmiss = []
for tid in sorted(frag):
    fname, frags = frag[tid]
    for fr in frags:
        hit = norm(fr).lower() in norm(files.get(fname, "")).lower()
        print("  %-4s %-22s %-70s %s" % (tid, fname, fr[:70], "found" if hit else "not-found"))
        if not hit: fragmiss.append((tid, fr))

print("\nSUMMARY: parse OK; structure problems=%d; anchored-quote misses=%d %s; grouped fragments not found=%d"
      % (len(problems), len(misses), misses, len(fragmiss)))
sys.exit(0)
