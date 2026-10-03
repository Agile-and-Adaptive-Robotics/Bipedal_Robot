"""What ARE the duplicated elements? For every duplicated ID, print each owning
element's tag/ClassName/Name so we can see which collisions are real damage
(e.g. a Link and a Gain sharing an ID) vs benign format repeats.

Run: D:/Anaconda/envs/myo/python.exe tools/probe_dups_20261002.py
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

APROJ = Path(__file__).resolve().parent.parent / "Walker_2_Layer_CPG_BilateralRG.aproj"
root = ET.parse(APROJ).getroot()

owners: dict[str, list[str]] = defaultdict(list)

def walk(e: ET.Element, path: str) -> None:
    eid = e.findtext("ID")
    if eid:
        cls = (e.findtext("ClassName") or "").strip().split(".")[-1]
        nm = (e.findtext("Name") or e.findtext("Text") or "").strip()
        owners[eid.strip()].append(f"{e.tag}/{cls} '{nm}'")
    for c in e:
        walk(c, path)

walk(root, "")

n_dup = 0
for eid, lst in sorted(owners.items()):
    if len(lst) > 1:
        n_dup += 1
        print(f"{eid[:13]}… x{len(lst)}: " + " | ".join(lst))
print(f"\ntotal duplicated IDs: {n_dup}")
