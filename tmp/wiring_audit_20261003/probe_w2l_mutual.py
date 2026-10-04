"""Audit probe 7: w2lvar IaIN/IBIN mutual directionality (same lazy
creation hazard as the stock full_rules branch)."""
import io
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, str(SPINAL))
os.chdir(SPINAL)
os.environ["AARL_NET"] = "w2lvar"

import muscle_map as MM  # noqa: E402
import build_network as bn  # noqa: E402
from build_network import ANTAGONIST  # noqa: E402

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")
net = bn.build(acts)
n = net.net
pop_names = [p["name"] for p in n.populations]
E = {(pop_names[c["source"]], pop_names[c["destination"]])
     for c in n.connections}
pairs = []
for act in acts:
    mi = MM.classify(act)
    for ant in ANTAGONIST.get(mi.groups[0], ()):
        for act2 in acts:
            mi2 = MM.classify(act2)
            if act2 != act and mi2.side == mi.side \
                    and mi2.groups[0] == ant:
                pairs.append((act, act2))
for fam in ("IaIN", "IBIN"):
    present = [(a, b) for a, b in pairs if (f"{fam}_{a}", f"{fam}_{b}") in E]
    hist = Counter()
    for a, b in pairs:
        k = tuple(sorted((a, b)))
        c = ((f"{fam}_{a}", f"{fam}_{b}") in E) \
            + ((f"{fam}_{b}", f"{fam}_{a}") in E)
        if c:
            hist[k] = c
    print(f"w2lvar {fam}: directed edges {len(present)} / {len(pairs)}; "
          f"unordered histogram {dict(Counter(hist.values()))}")
