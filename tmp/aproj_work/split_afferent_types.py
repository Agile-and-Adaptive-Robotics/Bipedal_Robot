"""Split the single 'Afferent HC Excite' synapse type in Walker_2_Layer_CPG_BilateralRG
into 20 uniquely named per-connection types (LH/RH symmetric instances share a type),
repoint all 39 instances, add the missing L heel->L Ank MN ext twin, delete the old type.
Text-level surgery: no reformatting of untouched XML. Physics identical (cloned type
params + untouched per-instance conductances)."""
import json
import re
import uuid
import xml.etree.ElementTree as ET

SRC = "W2LBRG.aproj"
DST = "W2LBRG_new.aproj"
OLD = "cafe0043-0000-4000-8000-900000000043"

def classify(src, tgt):
    s, t = src.lower(), tgt.lower()
    heel, toe = "heel contact" in s, "toe contact" in s
    aff = ("ia" in s.split()[-2:]) if False else None
    if heel:
        side_src = src[0]
        side_tgt = tgt[0]
        if "rg ext" in t: return "Heel_Ipsi_RGext_Excite"
        if "rg flx" in t: return "Heel_Contra_RGflx_Excite"
        if "hip pf ext" in t: return "Heel_Ipsi_HipPFext_Excite"
        if "knee pf ext" in t: return "Heel_Ipsi_KneePFext_Excite"
        if "hip mn ext" in t: return "Heel_Ipsi_HipMNext_Excite"
        if "knee mn ext" in t: return "Heel_Ipsi_KneeMNext_Excite"
        if "ank mn ext" in t: return "Heel_Ipsi_AnkMNext_Excite"
    if toe and "ank mn ext" in t: return "Toe_Ipsi_AnkMNext_Excite"
    if s.startswith(("l hip", "r hip", "l knee", "r knee")) or True:
        # muscle afferent sources: e.g. 'L Hip flx Ia' / 'R Knee ext Ib'
        parts = src.split(" ", 1)[1] if " " in src else src  # drop side letter
        aff = parts.split(" ")[-1]            # Ia / Ib / II
        grp = " ".join(parts.split(" ")[:-1])  # 'Hip flx' / 'Knee ext' ...
        tgt_short = tgt.split(" ", 1)[1] if " " in tgt else tgt  # 'RG flx' / 'Hip PF flx'
        grpu, tgu = grp.title().replace(" ", ""), tgt_short.replace(" ", "")
        return f"{aff}_{grpu}_to_{tgu}_Excite"
    return None

rows = json.load(open("aff_instances.json"))
raw = open(SRC, encoding="utf-8", errors="surrogateescape").read()

# --- locate + excise the old type block ---
idpos = raw.index(f"<ID>{OLD}</ID>")
tstart = raw.rindex("<Link>", 0, idpos)
tend = raw.index("</Link>", idpos) + len("</Link>")
old_block = raw[tstart:tend]
phys = [ln.strip() for ln in old_block.splitlines()
        if any(p in ln for p in ("EquilibriumPotential", "MaxSynapticConductance",
                                  "PreSynapticThreshold", "PreSynapticSaturationLevel"))]
print("old type physics params (cloned verbatim):")
for p in phys:
    print("  ", p)

# classify instances in document order, assign new type guids
new_guids = {}
order_guids = []
rows_by_pair = []
for r in rows:
    cls = classify(r["src"], r["tgt"])
    assert cls, f"unclassified: {r['src']} -> {r['tgt']}"
    if cls not in new_guids:
        n = len(new_guids) + 0x44
        new_guids[cls] = f"cafe00{n:02x}-0000-4000-8000-9000000000{n:02x}"
    order_guids.append(new_guids[cls])
    rows_by_pair.append((r["origin_guid"], r["dest_guid"], new_guids[cls], cls))
print(f"\nnew types ({len(new_guids)}):")
for k, v in new_guids.items():
    print(f"  {k}  {v}")

# --- build the replacement type blocks ---
def make_block(name, guid):
    b = old_block.replace(f"<ID>{OLD}</ID>", f"<ID>{guid}</ID>")
    b = b.replace("<Name>Afferent HC Excite</Name>", f"<Name>{name}</Name>")
    return b

type_blocks = "\r\n      ".join(make_block(k, v) for k, v in new_guids.items())

# excise old block, then repoint instances in document order
rest = raw[:tstart] + "@@TYPESLOT@@" + raw[tend:]
cnt = {"i": 0}
def sub(m):
    g = order_guids[cnt["i"]]
    cnt["i"] += 1
    return f"<SynapticTypeID>{g}</SynapticTypeID>"
rest = re.sub(r"<SynapticTypeID>" + re.escape(OLD) + r"</SynapticTypeID>", sub, rest)
assert cnt["i"] == 39, f"expected 39 instance repoints, got {cnt['i']}"
rest = rest.replace("@@TYPESLOT@@", type_blocks)
assert OLD not in rest, "old GUID still referenced"

# --- add the missing L heel -> L Ank MN ext twin (clone of the R twin) ---
r_twin = next(r for i, r in enumerate(rows) if i == 2)      # R heel -> R Ank MN ext
l_heel = next(r for i, r in enumerate(rows) if i == 29)["origin_guid"]   # L heel contact
l_ankmn = next(r for i, r in enumerate(rows) if i == 1)["dest_guid"]     # L Ank MN ext
ank_type = new_guids["Heel_Ipsi_AnkMNext_Excite"]

pos = 0
while True:
    ls = rest.find("<Link>", pos)
    if ls < 0:
        raise SystemExit("R twin block not found")
    le = rest.index("</Link>", ls) + len("</Link>")
    blk = rest[ls:le]
    if r_twin["origin_guid"] in blk and r_twin["dest_guid"] in blk:
        break
    pos = le
new_inst = blk
new_inst = re.sub(r"<ID>[^<]+</ID>", f"<ID>{str(uuid.uuid4())}</ID>", new_inst, count=1)
new_inst = new_inst.replace(f"<OriginID>{r_twin['origin_guid']}</OriginID>",
                            f"<OriginID>{l_heel}</OriginID>")
new_inst = new_inst.replace(f"<DestinationID>{r_twin['dest_guid']}</DestinationID>",
                            f"<DestinationID>{l_ankmn}</DestinationID>")
new_inst = re.sub(r"<SynapticTypeID>[^<]+</SynapticTypeID>",
                  f"<SynapticTypeID>{ank_type}</SynapticTypeID>", new_inst, count=1)
rest = rest[:le] + "\r\n      " + new_inst + rest[le:]

open(DST, "w", encoding="utf-8", errors="surrogateescape", newline="").write(rest)

# --- validate ---
root = ET.parse(DST).getroot()
def short(t):
    return t.split("}")[-1] if "}" in t else t
from collections import Counter
per = Counter()
total = 0
for el in root.iter():
    cm = {}
    for c in el:
        cm.setdefault(short(c.tag), c)
    tid = cm.get("SynapticTypeID")
    if tid is not None:
        total += 1
        per[(tid.text or "").strip()[:20]] += 1
print(f"\nVALIDATION: parse OK, total synapse instances = {total} (was 203, expect 204)")
old_refs = rest.count(OLD)
print(f"old-GUID references remaining: {old_refs}")
print(f"instances using new afferent types: {sum(v for k, v in per.items() if k.startswith('cafe004') or k.startswith('cafe005'))} (expect 40)")
