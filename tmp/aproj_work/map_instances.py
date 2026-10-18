"""Resolve synapse instances: source neuron, target (parent) neuron, type name,
per-instance params. Usage: map_instances.py <aproj> [type-filter]"""
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

DRAW_PROPS = {"ID", "Size", "Style", "Angle", "Filled", "AssemblyFile", "ClassName",
              "ZOrder", "Text", "BackMode", "DashStyle", "DrawColor", "DrawWidth",
              "Font", "ToolTip", "Url", "DestinationID", "OriginID", "AdjustDst",
              "AdjustOrg", "ArrowDestination", "ArrowMiddle", "ArrowOrigin", "Hidden",
              "Jump", "LineStyle", "OrthogonalDynamic", "OrientedText", "Selectable",
              "Stretchable", "X", "Y", "Z", "UserText", "Link", "Enabled", "Name"}

def short(tag):
    return tag.split("}")[-1] if "}" in tag else tag

def props(el):
    return {short(c.tag): (c.text or "").strip() for c in el}

def main(path, filt=None):
    root = ET.parse(path).getroot()
    parent = {c: p for p in root.iter() for c in p}
    by_guid = {}
    types = {}
    for el in root.iter():
        p = props(el)
        g = p.get("ID")
        if g:
            by_guid[g] = el
        if short(el.tag) == "SynapseTypes":
            for st in el:
                sp = props(st)
                if sp.get("ID"):
                    types[sp["ID"]] = sp.get("Name", "?")
                    by_guid[sp["ID"]] = st
    print(f"=== {path}  types={len(types)}")

    inst_params = Counter()
    per_type = defaultdict(list)
    for el in root.iter():
        p = props(el)
        tid = p.get("SynapticTypeID")
        if not tid:
            continue
        for k in p:
            if k not in ("SynapticTypeID",):
                inst_params[k] += 1
        # target: nearest ancestor with a Name that is a neuron-ish object
        tgt = el
        tgt_name = "?"
        chain = []
        cur = el
        while cur is not None and len(chain) < 4:
            cur = parent.get(cur)
            if cur is None:
                break
            cp = props(cur)
            nm = cp.get("Name")
            if nm:
                chain.append(nm)
                if len(chain) == 1:
                    tgt_name = nm
        src = p.get("SourceID", "?")
        src_el = by_guid.get(src)
        src_name = "?"
        if src_el is not None:
            sp = props(src_el)
            src_name = sp.get("Name", "?")
            pp = parent.get(src_el)
            if pp is not None:
                ppn = props(pp).get("Name")
                if ppn:
                    src_name = f"{ppn}/{src_name}"
        tname = types.get(tid, f"?{tid[:8]}")
        per_type[tname].append((src_name, tgt_name, p.get("SynapticConductance", "?")))
    print("instance property census:", dict(inst_params))
    print("\ninstances per type:")
    for t, lst in sorted(per_type.items(), key=lambda kv: -len(kv[1])):
        print(f"  {t}: {len(lst)}")
    if filt:
        print(f"\ninstances of types matching '{filt}':")
        for t, lst in per_type.items():
            if filt.lower() in t.lower():
                for s, tg, c in lst:
                    print(f"  [{t}] {s}  ->  {tg}   cond={c}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
