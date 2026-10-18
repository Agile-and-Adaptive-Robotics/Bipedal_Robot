"""Compare how heel/toe contact sensing is implemented in two AnimatLab models:
W2L BilateralRG vs Li. For each contact-labeled neuron: its class (NonSpiking vs
spiking), and any PhysicalToNodeAdapter feeding it (with the source body label)."""
import sys
import xml.etree.ElementTree as ET

def short(t):
    return t.split("}")[-1] if "}" in t else t

def cm(el):
    m = {}
    for c in el:
        m.setdefault(short(c.tag), c)
    return m

def label(el):
    c = cm(el)
    for k in ("Text", "UserText", "Name"):
        n = c.get(k)
        if n is not None and (n.text or "").strip():
            return (n.text or "").strip()
    return "?"

def analyze(path, needles):
    root = ET.parse(path).getroot()
    by_guid = {}
    adapters = []          # (origin_guid, dest_guid, class)
    for el in root.iter():
        c = cm(el)
        idel = c.get("ID")
        g = (idel.text or "").strip() if idel is not None and (idel.text or "").strip() else None
        cls = ""
        clsn = c.get("ClassName")
        if clsn is not None and clsn.text:
            cls = clsn.text.strip()
        if g:
            by_guid[g] = el
        if "Adapter" in cls or "ToNodeAdapter" in cls:
            org = c.get("OriginID"); dst = c.get("DestinationID")
            adapters.append(((org.text or "").strip() if org is not None else "?",
                             (dst.text or "").strip() if dst is not None else "?",
                             cls.split(".")[-1]))
    print(f"=== {path}")
    print(f"adapters found: {len(adapters)}")
    for needle in needles:
        for g, el in by_guid.items():
            lab = label(el)
            if needle.lower() in lab.lower():
                c = cm(el)
                cls = (c.get("ClassName").text or "?").split(".")[-1] if c.get("ClassName") is not None else "?"
                # find adapters feeding this node
                feeds = [(o, a) for o, d, a in adapters if d == g]
                src = "?"
                if feeds:
                    src = label(by_guid[feeds[0][0]]) if feeds[0][0] in by_guid else "?"
                print(f"  node '{lab}' class={cls} fed_by_adapter={len(feeds)} source='{src}'")

analyze(sys.argv[1] if len(sys.argv) > 1 else "W2LBRG.aproj", ["heel contact", "toe contact"])
print()
analyze("LI.aproj", ["contact", "heel", "touch"])
