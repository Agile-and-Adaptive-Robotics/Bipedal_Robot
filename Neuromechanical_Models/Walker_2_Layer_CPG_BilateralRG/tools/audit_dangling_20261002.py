"""Dangling-reference audit v2 (ElementTree) — 2026-10-02.

Checks every structural reference the AnimatLab save validator could trip on:
  1. every <Link> element's OriginID/DestinationID resolves to an existing element ID;
  2. every InLinks/OutLinks collection entry corresponds to a real <Link> element;
  3. every <Link> is registered in its origin's OutLinks AND destination's InLinks;
  4. every drawing <Tag> resolves to a Link or element ID;
  5. duplicate element IDs.

Run: D:/Anaconda/envs/myo/python.exe tools/audit_dangling_20261002.py
"""
from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

APROJ = Path(__file__).resolve().parent.parent / "Walker_2_Layer_CPG_BilateralRG.aproj"
tree = ET.parse(APROJ)
root = tree.getroot()

# --- pass 1: every element that directly owns an <ID> child -----------------
elem_info: dict[str, tuple[str, str, int]] = {}   # id -> (tag, classname, name)
dup: dict[str, int] = defaultdict(int)
all_elems: list[ET.Element] = []

def walk(e: ET.Element) -> None:
    all_elems.append(e)
    eid = e.findtext("ID")
    if eid:
        eid = eid.strip()
        dup[eid] += 1
        if eid not in elem_info:
            cls = (e.findtext("ClassName") or "").strip()
            name = (e.findtext("Name") or e.findtext("Text") or "").strip()
            elem_info[eid] = (e.tag, cls, name)
    for c in e:
        walk(c)

walk(root)

# --- links ------------------------------------------------------------------
links: list[tuple[str, str, str, str]] = []  # (id, origin, dest, classname)
for e in all_elems:
    if e.tag == "Link":
        cls = (e.findtext("ClassName") or "").strip()
        if "SynapseTypes" in cls:
            continue  # template definitions legitimately carry empty endpoints
        lid = (e.findtext("ID") or "").strip()
        org = (e.findtext("OriginID") or "").strip()
        dst = (e.findtext("DestinationID") or "").strip()
        links.append((lid, org, dst, cls))
link_by_id = {lid: (org, dst, cls) for lid, org, dst, cls in links if lid}

# --- collections ------------------------------------------------------------
coll: dict[str, list[tuple[str, str]]] = defaultdict(list)  # link id -> [(owner, kind)]
for e in all_elems:
    for kind in ("InLinks", "OutLinks"):
        blk = e.find(kind)
        if blk is None:
            continue
        owner = (e.findtext("ID") or "").strip() or e.tag
        for idel in blk.findall("ID"):
            coll[(idel.text or "").strip()].append((owner, kind))

problems: list[str] = []

for lid, org, dst, cls in links:
    label = f"link {lid or '(no ID)'} [{cls.split('.')[-1]}]"
    for role, ref in (("origin", org), ("destination", dst)):
        if not ref:
            problems.append(f"{label}: EMPTY {role}ID")
        elif ref not in elem_info:
            problems.append(f"{label}: {role} {ref} RESOLVES TO NOTHING")
    if lid:
        regs = coll.get(lid, [])
        kinds = {k for _, k in regs}
        if "In" not in kinds and "Out" not in kinds:
            problems.append(f"{label}: registered in NO InLinks/OutLinks collection")
        elif "In" not in kinds or "Out" not in kinds:
            problems.append(f"{label}: only in {sorted(kinds)} (expected In AND Out)")

for lid, regs in coll.items():
    if lid not in link_by_id:
        problems.append(f"collection entry {lid} (in {regs[0][0]} {regs[0][1]}) has NO <Link> element")

for e in all_elems:
    if e.tag == "Tag" or e.tag.endswith("Tag"):
        t = (e.text or "").strip()
        if t and t not in elem_info and t not in link_by_id:
            problems.append(f"drawing Tag {t} resolves to NOTHING")

for eid, n in dup.items():
    if n > 1:
        problems.append(f"DUPLICATE element ID {eid} defined {n} times")

print(f"elements-with-ID={len(elem_info)}  Link elements={len(links)}  collection entries={len(coll)}")
print(f"PROBLEMS: {len(problems)}")
for p in problems:
    print("  " + p)
sys.exit(0)
