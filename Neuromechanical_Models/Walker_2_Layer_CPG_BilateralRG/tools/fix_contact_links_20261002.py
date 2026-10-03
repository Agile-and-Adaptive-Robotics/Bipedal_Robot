"""Repair the unsavable Walker_2_Layer_CPG_BilateralRG.aproj (2026-10-02).

ROOT CAUSE (proven vs origin/AddingStepSensor_CoMorrow_stw reference):
the 9/16 contact surgery expressed the body->adapter binding as GUI
Links.Adapter elements ORIGINATING AT RigidBody contact bodies
(0221/0223/0225/0227). The loader cannot bind a behavioral link's origin to
a body, so a fresh load + save fails with "Link 'cafe0221-...' was missing
either an origin or destination node". The known-good reference never links
a body: the ADAPTER NODE itself carries <OriginID>(body) and
<DestinationID>(neuron) children.

FIX (text surgery only - the aproj contains CDATA drawings that XML writers
would mangle):
  1. backup to tools/backup_v5_prefix_20261002/
  2. insert <OriginID>/<DestinationID> into the four contact adapter nodes
     (before their <DataTypeID> line, matching the reference element order);
  3. delete the four body-origin GUI Links.Adapter elements;
  4. remove those four IDs from the adapter nodes' <InLinks>;
  5. leave all page drawings untouched (orphan drawing links are tolerated).
Asserts every precondition before touching anything.

Run: D:/Anaconda/envs/myo/python.exe tools/fix_contact_links_20261002.py
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
APROJ = HERE / "Walker_2_Layer_CPG_BilateralRG.aproj"
BACKUP = HERE / "tools" / "backup_v5_prefix_20261002"

text = APROJ.read_text(encoding="utf-8")
orig = text

def must(cond: bool, msg: str) -> None:
    if not cond:
        print(f"ABORT (no changes written): {msg}")
        sys.exit(1)

def find_block(t: str, start_re: str, end_tag: str) -> tuple[int, int]:
    m = re.search(start_re, t, re.M)
    must(m is not None, f"cannot locate {start_re[:60]}")
    s = m.start()
    e = t.index(end_tag, s) + len(end_tag)
    return s, e

# --- 1. resolve the four chains: body GUID -> illegal link -> adapter node -> downstream link -> neuron
SIDE = {
    "L heel": dict(body="foot_L_contact", adapter="L heel contact Adapter", neuron="L heel contact"),
    "L toe":  dict(body="toe_L_contact",  adapter="L toe contact Adapter",  neuron="L toe contact"),
    "R heel": dict(body="foot_R_contact", adapter="R heel contact Adapter", neuron="R heel contact"),
    "R toe":  dict(body="toe_R_contact",  adapter="R toe contact Adapter",  neuron="R toe contact"),
}

def guid_of_text(t: str, label: str, kind: str) -> str:
    """Resolve an element GUID by label via ElementTree (parse-only, no writing).
    kind: 'body' -> <RigidBody><Name>label; 'adapter' -> PhysicalToNodeAdapter with <Text>label;
          'neuron' -> ClassName .NonSpiking with <Text>label."""
    import xml.etree.ElementTree as ET
    root = ET.fromstring(t)
    for e in root.iter():
        eid = (e.findtext("ID") or "").strip()
        if not eid:
            continue
        cls = (e.findtext("ClassName") or "").strip()
        nm = (e.findtext("Name") or e.findtext("Text") or "").strip()
        if nm != label:
            continue
        if kind == "body" and e.tag == "RigidBody":
            return eid
        if kind == "adapter" and "PhysicalToNodeAdapter" in cls:
            return eid
        if kind == "neuron" and cls.endswith(".NonSpiking"):
            return eid
    must(False, f"cannot resolve GUID for '{label}' ({kind})")
    return ""

plan: dict[str, dict[str, str]] = {}
for side, spec in SIDE.items():
    body = guid_of_text(text, spec["body"], "body")
    adap = guid_of_text(text, spec["adapter"], "adapter")
    neuro = guid_of_text(text, spec["neuron"], "neuron")
    # find precisely by scanning all Adapter links for OriginID=body DestinationID=adapter
    ill_id = None
    for m in re.finditer(
        r"<ClassName>AnimatGUI\.DataObjects\.Behavior\.Links\.Adapter</ClassName>\s*\n"
        r"<ID>([0-9a-fA-F-]{36})</ID>(.*?)</Link>", text, re.S):
        body_block = m.group(2)
        o = re.search(r"<OriginID>([0-9a-fA-F-]{36})</OriginID>", body_block)
        d = re.search(r"<DestinationID>([0-9a-fA-F-]{36})</DestinationID>", body_block)
        if o and d and o.group(1) == body and d.group(1) == adap:
            ill_id = m.group(1)
            break
    must(ill_id is not None, f"{side}: no Adapter link {body[:8]}->{adap[:8]} found")
    # adapter node must already carry the correct body binding (reference pattern)
    am = re.search(rf"<ID>{adap}</ID>(.*?)</Node>", text, re.S)
    must(am is not None, f"{side}: adapter node block not found")
    blk = am.group(1)
    has_bind = f"<OriginID>{body}</OriginID>" in blk and f"<DestinationID>{neuro}</DestinationID>" in blk
    plan[side] = dict(body=body, adapter=adap, neuron=neuro, link=ill_id, bound=has_bind)
    print(f"{side}: body={body[:8]} adapter={adap[-4:]} neuron={neuro[-4:]} "
          f"illegal_link={ill_id[-4:]} node-binding={'present' if has_bind else 'MISSING'}")

# --- 2. backup
if not BACKUP.exists():
    BACKUP.mkdir(parents=True)
    shutil.copy2(APROJ, BACKUP / APROJ.name)
    print(f"backup written: {BACKUP / APROJ.name}")

# --- 3. apply edits (order: remove links, then InLinks entries, then insert node bindings)
for side, p in plan.items():
    # delete the illegal Link element (find its full <Link>...</Link> block)
    m = re.search(
        rf"<Link>\s*\n<AssemblyFile>AnimatGUI\.dll</AssemblyFile>\s*\n"
        rf"<ClassName>AnimatGUI\.DataObjects\.Behavior\.Links\.Adapter</ClassName>\s*\n"
        rf"<ID>{p['link']}</ID>.*?</Link>", text, re.S)
    must(m is not None, f"{side}: illegal link block vanished")
    text = text[:m.start()] + text[m.end():]

    # remove the link's ID from the adapter's InLinks
    inl = re.search(rf"(<InLinks>\s*\n)<ID>{p['link']}</ID>(\s*\n</InLinks>)", text)
    must(inl is not None, f"{side}: InLinks entry for {p['link'][-4:]} not found")
    text = text[:inl.start()] + "<InLinks>\n</InLinks>" + text[inl.end():]

    # ensure the adapter node carries the body binding (insert only if missing)
    if not p["bound"]:
        am = re.search(rf"<ID>{p['adapter']}</ID>(.*?)<DataTypeID>", text, re.S)
        must(am is not None, f"{side}: adapter insertion point not found")
        ins = f"<OriginID>{p['body']}</OriginID>\n<DestinationID>{p['neuron']}</DestinationID>\n"
        at = am.end(1)
        text = text[:at] + ins + text[at:]
        print(f"{side}: link removed, InLinks cleared, node binding INSERTED")
    else:
        print(f"{side}: link removed, InLinks cleared (node binding already present)")

# --- 4. sanity: well-formed, no trace of the four links outside drawings, count unchanged otherwise
import xml.etree.ElementTree as ET
from xml.dom import minidom
try:
    ET.fromstring(text)
except ET.ParseError as e:
    must(False, f"edited XML no longer well-formed: {e} - restoring backup")
for p in plan.values():
    n = text.count(p["link"])
    must(n == 1, f"link {p['link'][-4:]} still referenced {n} times (expected 1: drawing Tag only)")

APROJ.write_text(text, encoding="utf-8", newline="")
print(f"\nOK: {APROJ.name} repaired ({len(orig)} -> {len(text)} chars). "
      f"Reopen in AnimatLab2 and save - the four contact bindings now live on the adapter nodes "
      f"(reference pattern), and the body-origin GUI links are gone.")
