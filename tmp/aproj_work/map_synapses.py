"""Map AnimatLab synapse types and instances in an .aproj.

Structure discovered: synapse TYPES are elements under SynapseTypes (GUID in <ID> child);
synapse INSTANCES are children of the target neuron carrying <SynapticTypeID> and
<SynapticConductance>, with the source neuron referenced by <SourceID>.
"""
import sys
import xml.etree.ElementTree as ET

DRAW_PROPS = {"ID", "Size", "Style", "Angle", "Filled", "AssemblyFile", "ClassName",
              "ZOrder", "Text", "BackMode", "DashStyle", "DrawColor", "DrawWidth",
              "Font", "ToolTip", "Url", "DestinationID", "OriginID", "AdjustDst",
              "AdjustOrg", "ArrowDestination", "ArrowMiddle", "ArrowOrigin", "Hidden",
              "Jump", "LineStyle", "OrthogonalDynamic", "OrientedText", "Selectable",
              "Stretchable", "X", "Y", "Z", "UserText", "Link", "Enabled"}

def short(tag):
    return tag.split("}")[-1] if "}" in tag else tag

def props(el):
    return {short(c.tag): (c.text or "").strip() for c in el}

def build_maps(root):
    by_guid = {}
    synapse_types = []   # elements under SynapseTypes
    for el in root.iter():
        p = props(el)
        g = p.get("ID")
        if g:
            by_guid[g] = (short(el.tag), p.get("Name", "?"), el)
        if short(el.tag) == "SynapseTypes":
            for st in el:
                sp = props(st)
                if sp.get("ID"):
                    synapse_types.append(st)
                    by_guid[sp["ID"]] = (short(st.tag), sp.get("Name", "?"), st)
    return by_guid, synapse_types

def type_params(st):
    p = props(st)
    return {k: v for k, v in p.items()
            if k not in DRAW_PROPS and k not in ("Name",)}

def main(path, filter_type=None, list_all_types=False):
    root = ET.parse(path).getroot()
    by_guid, synapse_types = build_maps(root)
    print(f"=== {path}")
    print(f"synapse TYPES: {len(synapse_types)}")

    if list_all_types:
        for st in synapse_types:
            p = props(st)
            print(f"  TYPE {p.get('Name','?')}  guid={p.get('ID','?')[:8]}  {type_params(st)}")
        return

    tname = {}
    for st in synapse_types:
        p = props(st)
        tname[p.get("ID")] = p.get("Name", "?")

    # synapse instances: any element with a SynapticTypeID child
    n_inst = 0
    for el in root.iter():
        p = props(el)
        tid = p.get("SynapticTypeID")
        if not tid:
            continue
        n_inst += 1
        t = tname.get(tid, f"?{tid[:8]}")
        if filter_type and filter_type.lower() not in t.lower():
            continue
        # parent chain: the synapse element's parent is the target neuron
        src = p.get("SourceID", "?")
        stag, sname, _ = by_guid.get(src, ("", src[:8] if src != "?" else "?", None))
        # find enclosing neuron name: walk up via parent map
        print(f"  INST type='{t}' src='{sname}' cond={p.get('SynapticConductance','?')}"
              f" g={p.get('Gmax') or p.get('Gain') or ''} el={short(el.tag)}"
              f" parent={p.get('Name','?')}")
    print(f"synapse INSTANCES total: {n_inst}")

if __name__ == "__main__":
    a = sys.argv
    if len(a) > 2 and a[2] == "--types":
        main(a[1], list_all_types=True)
    else:
        main(a[1], a[2] if len(a) > 2 else None)
