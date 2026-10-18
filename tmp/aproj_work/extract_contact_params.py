"""Extract from LI.aproj: the ground-contact bodies (parent, position, size, contact
flags) and the adapters feeding the contact spiking neurons (type, gains).
From W2LBRG.aproj: foot-body names (attachment parents) + the 4 contact neurons'
diagram pages and Location coordinates."""
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

def guid(el):
    c = cm(el)
    idel = c.get("ID")
    return (idel.text or "").strip() if idel is not None and idel.text else None

def interesting_props(el):
    out = {}
    for c in el:
        t = short(c.tag)
        if t in ("Location", "Size", "DrawColor", "Font", "ToolTip", "Url", "Text",
                 "UserText", "Alignment", "AutoSize", "BackMode", "DashStyle",
                 "ShadowStyle", "ZOrder", "InLinkable", "OutLinkable", "LabelEdit",
                 "Grayscale", "Editable", "Selectable", "Resizable", "Movable",
                 "ArrowDestination", "ArrowMiddle", "ArrowOrigin", "Jump", "LineStyle",
                 "OrthogonalDynamic", "OrientedText", "AdjustDst", "AdjustOrg",
                 "Hidden", "Stretchable", "PortUnion", "From", "Port", "Bottom", "Top"):
            continue
        v = (c.text or "").strip()
        if v or c.attrib:
            out[t] = v if v else c.attrib
    return out

def dump(path, label_needles, body_needles):
    root = ET.parse(path).getroot()
    by_guid = {}
    parent = {c: p for p in root.iter() for c in p}
    for el in root.iter():
        g = guid(el)
        if g:
            by_guid[g] = el
    print(f"=== {path}")
    # adapters
    for el in root.iter():
        cls = (cm(el).get("ClassName").text or "") if cm(el).get("ClassName") is not None else ""
        if "Adapter" in cls:
            c = cm(el)
            dst = (c.get("DestinationID").text or "").strip() if c.get("DestinationID") is not None else ""
            tgt = label(by_guid[dst]) if dst in by_guid else "?"
            if any(n.lower() in tgt.lower() for n in label_needles):
                org = (c.get("OriginID").text or "").strip() if c.get("OriginID") is not None else ""
                src = label(by_guid[org]) if org in by_guid else "?"
                print(f"ADAPTER {cls.split('.')[-1]} '{label(el)}' src='{src}' -> '{tgt}' props={interesting_props(el)}")
    # bodies
    for el in root.iter():
        lab = label(el)
        if any(n.lower() in lab.lower() for n in body_needles):
            cls = (cm(el).get("ClassName").text or "?") if cm(el).get("ClassName") is not None else short(el.tag)
            par = parent.get(el)
            plab = label(par) if par is not None else "?"
            print(f"BODY '{lab}' class={cls.split('.')[-1]} parent='{plab}' props={interesting_props(el)}")
    return by_guid, parent

dump("LI.aproj", ["contact"], ["foot ground", "toe ground"])
print()
by_guid, parent = dump("W2LBRG.aproj", ["heel contact", "toe contact"], ["foot", "ankle", "shank", "tibia"])
# contact neurons' page + location
root2 = ET.parse("W2LBRG.aproj").getroot()
parent2 = {c: p for p in root2.iter() for c in p}
for el in root2.iter():
    lab = label(el)
    if lab in ("L heel contact", "R heel contact", "L toe contact", "R toe contact"):
        loc = cm(el).get("Location")
        par = parent2.get(el)
        chain = []
        p = par
        while p is not None and len(chain) < 3:
            l = label(p)
            if l != "?":
                chain.append(l)
            p = parent2.get(p)
        print(f"NODE '{lab}' loc={loc.attrib if loc is not None else '?'} parents={chain}")
