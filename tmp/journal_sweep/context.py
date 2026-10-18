# -*- coding: utf-8 -*-
"""Cross-check flagged lines vs baseline; print full context for judgment items."""
import re, os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
ROOT = r"D:\Github\Bipedal_Robot"
CH = os.path.join(ROOT, r"Documentation\Reports and Papers\Dissertation\ProofFinal\chapters")
BASE = os.path.join(ROOT, r"tmp\diss_baseline_20261017")
def load(p):
    for enc in ("utf-8","latin-1"):
        try: return open(p, encoding=enc).read()
        except UnicodeDecodeError: continue
def squash(s): return re.sub(r'\s+','',s)

items = [
 ("15-background", [8,66,111]),
 ("20-muscles", [74,94,328,369]),
 ("30-torque", [18,32]),
 ("60-sns", [285,316,333,366]),
 ("93-AppendixB", [87,110,171]),
 ("94-AppendixC", [34,91,95,126,138,139,140,157,207,209,236,252,114,121]),
 ("70-overall", [42]),
 ("20-muscles", [295,310,410]),
 ("60-sns", [135]),
]
from collections import defaultdict
byfile = defaultdict(list)
for f, ls in items: byfile[f].extend(ls)
for f, ls in byfile.items():
    cur = load(os.path.join(CH, f+".tex")); base = load(os.path.join(BASE, f+".tex"))
    sqb = squash(base)
    curl = cur.splitlines()
    print("="*100); print("FILE", f)
    for n in sorted(set(ls)):
        ln = curl[n-1] if n<=len(curl) else "<line beyond EOF>"
        changed = squash(ln) not in sqb
        print(f"\n--- L{n} [{'CHANGED/NEW vs baseline' if changed else 'in baseline (grandfathered?)'}]")
        print(ln)
print("="*100)
# diff-added line listing for 30-torque around L32, and a compact added-lines census per chapter
print("\n\nADDED-LINE CENSUS (added lines in unified diff, per chapter):")
CHAPTERS = ["10-introduction","15-background","20-muscles","30-torque","50-database",
            "60-sns","70-overall","80-conclusion","93-AppendixB","94-AppendixC"]
for name in CHAPTERS:
    cur = load(os.path.join(CH, name+".tex")); base = load(os.path.join(BASE, name+".tex"))
    from collections import Counter
    cb, bb = Counter(squash(l) for l in base.splitlines()), Counter()
    added = []
    cb = Counter(squash(l) for l in base.splitlines())
    for l in cur.splitlines():
        s = squash(l)
        if s and cb[s]>0: cb[s]-=1
        else: added.append(l)
    print(f"  {name}: {len(added)} added/changed lines")
