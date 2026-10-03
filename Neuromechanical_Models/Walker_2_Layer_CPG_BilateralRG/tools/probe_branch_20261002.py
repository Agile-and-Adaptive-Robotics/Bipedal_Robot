"""How does the KNOWN-GOOD branch W2L bind contact bodies to adapters?
Extract from /tmp/w2l_branch_ref.aproj:
  - every Links.Adapter whose origin resolves to a RigidBody (vs Node)
  - the first contact adapter node's full config (Source binding children)
Run: D:/Anaconda/envs/myo/python.exe tools/probe_branch_20261002.py
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

REF = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\w2l_branch_ref.aproj")
if not REF.exists():
    REF = Path("/tmp/w2l_branch_ref.aproj")
root = ET.parse(REF).getroot()

real: dict[str, ET.Element] = {}
def idx(e: ET.Element) -> None:
    if e.tag != "ID" and e.find("ID") is not None:
        eid = (e.findtext("ID") or "").strip()
        real.setdefault(eid, e)
    for c in e:
        idx(c)
idx(root)

def kind(e: ET.Element) -> str:
    if e is None:
        return "UNRESOLVED"
    cls = (e.findtext("ClassName") or "").strip().split(".")[-1]
    nm = (e.findtext("Name") or e.findtext("Text") or "").strip()
    return f"{e.tag}/{cls} '{nm}'"

hist: Counter[str] = Counter()
body_links: list[str] = []
for e in root.iter("Link"):
    cls = (e.findtext("ClassName") or "").strip()
    if "Links.Adapter" not in cls:
        continue
    org = (e.findtext("OriginID") or "").strip()
    dst = (e.findtext("DestinationID") or "").strip()
    oel, del_ = real.get(org), real.get(dst)
    ok, dk = kind(oel), kind(del_)
    hist[f"origin {'RIGIDBODY' if ok.startswith('RigidBody') else 'node/other'} -> dest {'RIGIDBODY' if dk.startswith('RigidBody') else 'node/other'}"] += 1
    if ok.startswith("RigidBody") or dk.startswith("RigidBody"):
        body_links.append(f"  { (e.findtext('ID') or '?')[-6:] }: {ok} -> {dk}")

print("== Adapter link origin/dest kinds (branch reference) ==")
for k, v in hist.most_common():
    print(f"  {v:3d} x  {k}")
print("\n== body-touching adapter links ==")
print("\n".join(body_links) or "  NONE")

# first contact adapter node config
print("\n== first PhysicalToNodeAdapter node (full child list) ==")
for e in root.iter():
    if "PhysicalToNodeAdapter" in (e.findtext("ClassName") or ""):
        nm = e.findtext("Text") or e.findtext("Name") or ""
        kids = [c.tag for c in e]
        print(f"'{nm}' children: {kids}")
        for c in e:
            if c.tag in ("SourceID", "DataTypeID", "TargetDataTypeID", "InLinks", "OutLinks", "Gain", "BodyID", "OriginID", "DestinationID"):
                txt = (c.text or "").strip()
                sub = [g.tag for g in c]
                print(f"   {c.tag}: {txt or sub}")
        break
