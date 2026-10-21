"""Antagonist-pair coverage check (2026-10-20 audit): for every pair in
build_network.ANTAGONIST, count the IaIN->MN and IIIN->MN edges syn6
builds, on the full 92-actuator set. One line per pair; the sum must
equal the 436+436 the edge dump reports."""
from __future__ import annotations

import io
import os
import sys
from collections import defaultdict

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
spinal = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, spinal)
os.chdir(spinal)

import muscle_map as MM                  # noqa: E402
import build_network as BN                # noqa: E402

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")
muscles = {a: MM.classify(a) for a in acts}

pair_counts = defaultdict(int)
muscle_count = defaultdict(int)
for act, mi in muscles.items():
    muscle_count[mi.groups[0]] += 1
for act, mi in muscles.items():
    for ant in BN.ANTAGONIST.get(mi.groups[0], ()):
        for act2, mi2 in muscles.items():
            if mi2.side == mi.side and mi2.groups[0] == ant:
                pair_counts[(mi.groups[0], ant)] += 1

tot = 0
print(f"{'pair':<24} {'n_src':>6} {'n_dst':>6} {'ordered edges':>14}")
for (g, ant) in sorted(pair_counts):
    n = pair_counts[(g, ant)]
    tot += n
    print(f"{g:>10} -> {ant:<12} {muscle_count[g]:>6} "
          f"{muscle_count[ant]:>6} {n:>14}")
print(f"total ordered antagonist pairs = {tot} "
      f"(IaIN->MN {tot}, IIIN->MN {tot}, IaIN<->IaIN {tot}, "
      f"IBIN<->IBIN {tot})")
