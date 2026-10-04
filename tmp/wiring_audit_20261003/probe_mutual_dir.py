"""Audit probe 6: quantify IaIN<->IaIN and IBIN<->IBIN mutual-inhibition
directionality in the STOCK build under full_rules (lazy IN creation
ordering hazard), and confirm muscle wiring order."""
import io
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, str(SPINAL))
os.chdir(SPINAL)
assert os.environ.get("AARL_NET", "") == ""

import muscle_map as MM  # noqa: E402
import params as P  # noqa: E402
import build_network as bn  # noqa: E402
from build_network import ANTAGONIST  # noqa: E402

P.G["full_rules"] = 1.0
P.G["ia_in"] = 0.825

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")
print("wiring order (first 12):", acts[:12])
net = bn.build(acts)
n = net.net
pop_names = [p["name"] for p in n.populations]
E = set()
for c in n.connections:
    E.add((pop_names[c["source"]], pop_names[c["destination"]]))

# expected ordered antagonist pairs
pairs = []
for act in acts:
    mi = MM.classify(act)
    for ant in ANTAGONIST.get(mi.groups[0], ()):
        for act2 in acts:
            mi2 = MM.classify(act2)
            if act2 != act and mi2.side == mi.side \
                    and mi2.groups[0] == ant:
                pairs.append((act, act2))
print(f"ordered antagonist muscle pairs: {len(pairs)}")

iain_fwd = [(a, b) for a, b in pairs if (f"IaIN_{a}", f"IaIN_{b}") in E]
ibin_fwd = [(a, b) for a, b in pairs if (f"IBIN_{a}", f"IBIN_{b}") in E]
print(f"IaIN directed mutual edges present: {len(iain_fwd)} / {len(pairs)}")
print(f"IBIN directed mutual edges present: {len(ibin_fwd)} / {len(pairs)}")
# per unordered pair: count directions present
from collections import Counter  # noqa: E402
ui = Counter()
for a, b in pairs:
    key = tuple(sorted((a, b)))
    if (f"IaIN_{a}", f"IaIN_{b}") in E:
        ui[key] += 1
    if (f"IaIN_{b}", f"IaIN_{a}") in E:
        ui[key] += 1
hist = Counter(ui.values())
print(f"IaIN unordered-pair direction histogram (1=one-way, 2=mutual): "
      f"{dict(hist)}")
ub = Counter()
for a, b in pairs:
    key = tuple(sorted((a, b)))
    if (f"IBIN_{a}", f"IBIN_{b}") in E:
        ub[key] += 1
    if (f"IBIN_{b}", f"IBIN_{a}") in E:
        ub[key] += 1
print(f"IBIN unordered-pair direction histogram: {dict(Counter(ub.values()))}")
# which direction exists: later-wired -> earlier-wired?
order = {a: i for i, a in enumerate(acts)}
late2early = sum(1 for (s, d) in iain_fwd if order[s] > order[d])
print(f"IaIN edges pointing later->earlier in wiring order: {late2early}"
      f" / {len(iain_fwd)}")
late2early_b = sum(1 for (s, d) in ibin_fwd if order[s] > order[d])
print(f"IBIN edges pointing later->earlier in wiring order: "
      f"{late2early_b} / {len(ibin_fwd)}")
missing = [(a, b) for a, b in pairs
           if (f"IaIN_{a}", f"IaIN_{b}") not in E][:10]
print("sample missing IaIN directions:", missing)
