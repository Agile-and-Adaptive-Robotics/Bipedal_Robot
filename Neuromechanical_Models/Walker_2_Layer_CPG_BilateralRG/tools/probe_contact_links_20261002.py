"""Focused probe: the four contact Adapter links and their endpoint registration.

For each GUI Links.Adapter from a contact body -> adapter node, report:
  - origin/dest element type + name
  - does the DEST node's InLinks list the link?
  - does the ORIGIN body's OutLinks list the link (or does it even have OutLinks)?
  - the parent chain where the Link element sits

Run: D:/Anaconda/envs/myo/python.exe tools/probe_contact_links_20261002.py
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

APROJ = Path(__file__).resolve().parent.parent / "Walker_2_Layer_CPG_BilateralRG.aproj"
root = ET.parse(APROJ).getroot()

parent: dict[int, ET.Element] = {}
def build_parents(e: ET.Element) -> None:
    for c in e:
        parent[id(c)] = e
        build_parents(c)
build_parents(root)

def chain(e: ET.Element, depth: int = 4) -> str:
    out = []
    cur = e
    for _ in range(depth):
        p = parent.get(id(cur))
        if p is None:
            break
        label = p.tag
        pid = p.findtext("ID") or ""
        pname = p.findtext("Name") or p.findtext("Text") or ""
        if pid or pname:
            label += f"[{pname or pid[:8]}]"
        out.append(label)
        cur = p
    return " < ".join(out)

by_id: dict[str, ET.Element] = {}
def index(e: ET.Element) -> None:
    eid = e.findtext("ID")
    if eid:
        by_id.setdefault(eid.strip(), e)
    for c in e:
        index(c)
index(root)

TARGETS = [
    "cafe0221-0000-4000-8000-900000000221",
    "cafe0223-0000-4000-8000-900000000223",
    "cafe0225-0000-4000-8000-900000000225",
    "cafe0227-0000-4000-8000-900000000227",
]

for t in TARGETS:
    lk = by_id.get(t)
    print("=" * 70)
    if lk is None:
        print(f"{t}: NO LINK ELEMENT")
        continue
    cls = (lk.findtext("ClassName") or "").strip()
    org = (lk.findtext("OriginID") or "").strip()
    dst = (lk.findtext("DestinationID") or "").strip()
    print(f"link {t[-4:]} [{cls.split('.')[-1]}]  origin={org[-8:]} dest={dst[-8:]}")
    print(f"  parent chain: {chain(lk)}")
    for role, ref in (("ORIGIN", org), ("DEST", dst)):
        el = by_id.get(ref)
        if el is None:
            print(f"  {role}: {ref} NOT FOUND")
            continue
        cls2 = (el.findtext("ClassName") or el.tag).split(".")[-1]
        nm = el.findtext("Name") or el.findtext("Text") or ""
        inl = el.find("InLinks")
        outl = el.find("OutLinks")
        in_ids = [c.text for c in inl.findall("ID")] if inl is not None else None
        out_ids = [c.text for c in outl.findall("ID")] if outl is not None else None
        print(f"  {role}: {cls2} '{nm}' ({ref[:8]})")
        print(f"    InLinks : {'ABSENT' if in_ids is None else [i[-4:] if i else '?' for i in in_ids]}")
        print(f"    OutLinks: {'ABSENT' if out_ids is None else [i[-4:] if i else '?' for i in out_ids]}")

# also dump the four contact bodies' child tags for completeness
print("=" * 70)
for name in ("foot_L_contact", "toe_L_contact", "foot_R_contact", "toe_R_contact"):
    hits = [e for e in by_id.values() if (e.findtext("Name") or "").strip() == name]
    for el in hits:
        kids = sorted({c.tag for c in el})
        print(f"{name}: children={kids}")
