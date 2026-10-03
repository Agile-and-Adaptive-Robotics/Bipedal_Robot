"""Parent-container map for every Links.Adapter in the BilateralRG aproj:
which parent element each GUI Adapter link sits under, and where its origin
points (body vs node). Contact links 0221/0223/0225/0227 flagged.

Run: D:/Anaconda/envs/myo/python.exe tools/probe_parents_20261002.py
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

APROJ = Path(__file__).resolve().parent.parent / "Walker_2_Layer_CPG_BilateralRG.aproj"
root = ET.parse(APROJ).getroot()

parent: dict[int, ET.Element] = {}
def build(e: ET.Element) -> None:
    for c in e:
        parent[id(c)] = e
        build(c)
build(root)

def plabel(e: ET.Element) -> str:
    nm = e.findtext("Name") or e.findtext("Text") or e.findtext("ClassName") or ""
    nm = nm.split(".")[-1].strip()
    return f"{e.tag}:{nm}"

# element-with-ID index (exclude bare <ID> entries and collection containers)
real: dict[str, ET.Element] = {}
def idx(e: ET.Element) -> None:
    if e.tag not in ("ID",) and e.find("ID") is not None:
        eid = (e.findtext("ID") or "").strip()
        real.setdefault(eid, e)
    for c in e:
        idx(c)
idx(root)

CONTACT = {
    "cafe0221-0000-4000-8000-900000000221",
    "cafe0223-0000-4000-8000-900000000223",
    "cafe0225-0000-4000-8000-900000000225",
    "cafe0227-0000-4000-8000-900000000227",
}

counts: Counter[str] = Counter()
for e in root.iter("Link"):
    cls = (e.findtext("ClassName") or "").strip()
    if "Links.Adapter" not in cls:
        continue
    lid = (e.findtext("ID") or "").strip()
    org = (e.findtext("OriginID") or "").strip()
    p = parent.get(id(e))
    chain = plabel(p) + " < " + plabel(parent.get(id(p))) if p is not None and parent.get(id(p)) is not None else (plabel(p) if p is not None else "?")
    oel = real.get(org)
    otype = plabel(oel) if oel is not None else "UNRESOLVED"
    tag = "CONTACT " if lid in CONTACT else ""
    counts[chain.split(" < ")[0] + " || " + otype] += 1
    if lid in CONTACT:
        print(f"{tag}{lid[-4:]}: parent={chain}  origin({org[:8]})={otype}")

print("\n== parent || origin-type histogram (all Adapter links) ==")
for k, v in counts.most_common():
    print(f"  {v:3d} x  {k}")
