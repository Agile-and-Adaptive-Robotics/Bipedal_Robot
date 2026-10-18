"""Analyze AnimatLab .aproj synapse types: dump afferent-named synapses with
source/target resolution, and the K&A HF2l-3 reference params from Biped_2xCPG_wSubs."""
import sys
import xml.etree.ElementTree as ET

def load(path):
    tree = ET.parse(path)
    return tree.getroot()

def index_organisms(root):
    """Map every element with an ID attribute to (tag, name)."""
    idx = {}
    for el in root.iter():
        eid = el.attrib.get("ID")
        if eid is not None:
            idx[eid] = (el.tag, el.attrib.get("Name", "?"))
    return idx

def short(tag):
    return tag.split("}")[-1] if "}" in tag else tag

def dump_synapse(el, idx, indent="  "):
    src = el.attrib.get("SourceID") or el.attrib.get("FromID") or "?"
    dst = el.attrib.get("TargetID") or el.attrib.get("ToID") or "?"
    stag, sname = idx.get(src, (src, "?"))
    ttag, tname = idx.get(dst, (dst, "?"))
    params = {c.attrib.get("Name", short(c.tag)): (c.text or "").strip()
              for c in el if not c.attrib.get("ID")}
    print(f"{indent}{el.attrib.get('Name','?')} [{short(el.tag)}]  {sname} -> {tname}")
    if params:
        print(f"{indent}  params: {params}")
    for c in el:
        if c.attrib.get("ID"):
            print(f"{indent}  child[{short(c.tag)}] Name={c.attrib.get('Name')}")

def main(path, needle, show_children=False):
    root = load(path)
    idx = index_organisms(root)
    print(f"=== {path}: elements={len(idx)} ===")
    # every element whose tag mentions Synapse or Stimulus
    syn_tags = {}
    for el in root.iter():
        t = short(el.tag)
        if "Synapse" in t:
            syn_tags[t] = syn_tags.get(t, 0) + 1
    print("synapse tags:", syn_tags)
    for el in root.iter():
        t = short(el.tag)
        if "Synapse" not in t:
            continue
        name = el.attrib.get("Name", "")
        if needle.lower() in name.lower():
            dump_synapse(el, idx)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
