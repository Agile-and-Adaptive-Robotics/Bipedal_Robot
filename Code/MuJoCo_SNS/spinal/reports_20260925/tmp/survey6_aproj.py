# -*- coding: utf-8 -*-
"""Survey aproj: full tree with types, boxes with size/mass, joints, muscles."""
import xml.etree.ElementTree as ET

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()

def txt(el, tag):
    c = el.find(tag)
    return (c.text or "").strip() if c is not None else None

def act(el, tag, attr="Actual"):
    if el is None:
        return None
    c = el.find(tag)
    return c.get(attr) if c is not None else None

org = root.find(".//Organism")

def dump(rb, depth=0):
    pad = "  " * depth
    nm = txt(rb, "Name"); t = txt(rb, "Type")
    cid = (txt(rb, "ID") or "")[:8]
    extra = ""
    if t == "Box":
        sz = rb.find("Size")
        extra = f" size=({act(sz,'Length')},{act(sz,'Width')},{act(sz,'Height')}) mass={act(rb,'Mass')}g freeze={txt(rb,'Freeze')} matID={txt(rb,'MaterialTypeID')}"
        lp = rb.find("LocalPosition")
        if lp is not None:
            extra += f" pos=({act(lp,'X')},{act(lp,'Y')},{act(lp,'Z')})"
        ro = rb.find("Rotation")
        if ro is not None:
            extra += f" rot=({act(ro,'X')},{act(ro,'Y')},{act(ro,'Z')})"
        com = rb.find("COM")
        if com is not None:
            extra += f" com=({act(com,'X')},{act(com,'Y')},{act(com,'Z')})"
    elif t and "Muscle" in t:
        at = [a.text for a in rb.findall("Attachments/AttachID")]
        extra = f" maxT={act(rb,'MaximumTension')} Kse={act(rb,'Kse')} Kpe={act(rb,'Kpe')} B={act(rb,'B')} applyT={txt(rb,'ApplyTension')}"
        lt = rb.find("LengthTension")
        if lt is not None:
            extra += f" restLen={act(lt,'RestingLength')} Lwidth={act(lt,'Lwidth')}"
        extra += f" attach={at}"
        lp = rb.find("LocalPosition")
        if lp is not None:
            extra += f" pos=({act(lp,'X')},{act(lp,'Y')},{act(lp,'Z')})"
    elif t == "Attachment":
        lp = rb.find("LocalPosition")
        extra = f" pos=({act(lp,'X')},{act(lp,'Y')},{act(lp,'Z')})" if lp is not None else " (no pos)"
    elif t and "StretchReceptor" in t:
        extra = f" applyT={txt(rb,'ApplyTension')} maxT={act(rb,'MaximumTension')}"
    elif t == "Spring":
        extra = f" len={act(rb,'Length')}"
    # joints are direct children Joint of RigidBody
    for jt in rb.findall("Joint"):
        lp = jt.find("LocalPosition"); ro = jt.find("Rotation")
        print(f"{pad}   >JOINT {txt(jt,'Name')} pos=({act(lp,'X')},{act(lp,'Y')},{act(lp,'Z')}) rot=({act(ro,'X')},{act(ro,'Y')},{act(ro,'Z')})")
        el = jt.find("EnableLimits")
        if el is not None:
            print(f"{pad}      EnableLimits={el.text} Lower={act(jt,'LowerLimit')} Upper={act(jt,'UpperLimit')}")
    print(f"{pad}{nm} [{t}] id={cid}{extra}")
    cb = rb.find("ChildBodies")
    if cb is not None:
        for ch in cb:
            dump(ch, depth + 1)

for ch in org:
    if ch.tag == "RigidBody":
        dump(ch)
    elif ch.tag == "Structure":
        print("STRUCTURE:", txt(ch, "Name"))
        for c in ch:
            if c.tag == "RigidBody":
                dump(c, 1)
