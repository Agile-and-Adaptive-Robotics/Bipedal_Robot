import json
import os
import re

RC = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C"

cites = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_cites.json", encoding="utf-8"))
gid = "recFA1PeRSPTv9POu"  # Grillner 1985
out = cites.get(gid, [])
citers = [k for k, v in cites.items() if gid in v]
print("Grillner1985 cited-by (in-corpus):", len(citers), "| its reference list length:", len(out))

# pulley surface figure files
for root, dirs, files in os.walk(os.path.join(RC, "figs")):
    for f in files:
        if "ulley" in f or "Margin" in f or "Surface" in f:
            print("pulley fig:", os.path.relpath(os.path.join(root, f), RC))

# case-insensitive author check
for f in sorted(os.listdir(os.path.join(RC, "chapters"))):
    if f.endswith(".tex"):
        s = open(os.path.join(RC, "chapters", f), encoding="utf-8", errors="replace").read()
        for m in re.finditer(r"(?i)the author", s):
            i = m.start()
            print("PERSON:", f, "::", s[max(0, i - 45) : i + 45].replace("\n", " "))

# numeric ledger spot checks
nc = json.load(open(os.path.join(RC, "review", "numeric_change_review.json"), encoding="utf-8"))
print("numeric_change_review entries:", len(nc) if hasattr(nc, "__len__") else "?")
blob = json.dumps(nc)
for probe in ("4.58", "0.89", "8816", "8761", "0.525", "50.6", "28.8", "259.9", "1,393", "1393", "55.6"):
    print(f"ledger mentions {probe!r}:", probe in blob)
