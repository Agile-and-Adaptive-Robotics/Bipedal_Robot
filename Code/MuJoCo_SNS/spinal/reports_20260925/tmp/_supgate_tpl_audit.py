"""SUPERVISOR GATE: independent audit of connectome_templates.json
bilateralrg entry (the four claimed generation fixes) + template counts.
Read-only.
"""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
P = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json"
d = json.load(open(P, encoding="utf-8"))
print("template keys:", sorted(d.keys()))
b = d["bilateralrg"]
nodes, edges = b["nodes"], b["edges"]
print(f"bilateralrg: {len(nodes)} nodes / {len(edges)} edges "
      f"(claims: 109/222; pre-fix 105/226)")


def lbl(n):
    return n.get("label", "")


# error 1: V3 must target the contralateral RG ext IN, NOT the RG ext HC
v3e = [e for e in edges if str(e.get("from", "")).startswith("V3_")
       or str(e.get("from", "")) in ("L V3", "R V3")]
print("V3 out-edges:", [(e["from"], e["to"], e.get("sign"),
                         e.get("tag")) for e in v3e])
old = [e for e in edges if e["tag"] == "v3_SynAmp0.1_weak"]
print("old wrong tag v3_SynAmp0.1_weak present:", len(old))
new = [e for e in edges if e["tag"] == "v3_to_contra_InE"]
print("v3_to_contra_InE edges:", len(new))

# error 2: contact must NOT hit ankle MN ext
contact = [e for e in edges if e.get("tag") == "contact_C20"]
ank = [e for e in contact if "Ank" in str(e["to"])]
print(f"contact_C20 edges: {len(contact)} (claims 20); ankle targets: "
      f"{[(e['from'], e['to']) for e in ank]}")

# error 3: afferent-HC excite edges: 24, own-joint, incl. RG HC targets
aff = [e for e in edges if e.get("tag") == "aff_HC_excite"]
print(f"aff_HC_excite edges: {len(aff)} (claims 24)")
sides = {}
for e in aff:
    sides.setdefault(e["from"][0], set()).add(str(e["to"]))
for s in sorted(sides):
    print(f"  side {s} targets: {sorted(sides[s])}")
sn2 = [n for n in nodes if "II" in lbl(n) and lbl(n).startswith("S J")]
print("SN-II relay nodes:", sorted(lbl(n) for n in sn2))

# error 4: R RG<->RG direct pair present
rg = [e for e in edges if "RG" in str(e.get("from", ""))
      and "RG" in str(e.get("to", ""))]
print("all RG->RG edges:",
      [(e["from"], e["to"], e.get("sign"), e.get("tag")) for e in rg])
ii = [lbl(n) for n in nodes if " II" in lbl(n) or lbl(n).endswith("II")]
print("II-containing node labels:", sorted(set(ii)))
sn = [lbl(n) for n in nodes if str(lbl(n)).startswith("S ")]
print("'S '-prefixed (afferent) labels:", sorted(set(sn)))
