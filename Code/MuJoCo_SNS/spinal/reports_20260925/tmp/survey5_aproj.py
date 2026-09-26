# -*- coding: utf-8 -*-
"""Survey aproj: hierarchical RigidBody tree + one full muscle dump."""
import xml.etree.ElementTree as ET

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()

def txt(el, tag):
    c = el.find(tag)
    return (c.text or "").strip() if c is not None else None

def trip(el, tag):
    c = el.find(tag)
    if c is None:
        return None
    return (c.get("Value"), c.get("Scale"), c.get("Actual"))

org = root.find(".//Organism")
body = org.find("RigidBody")

def dump_rb(rb, depth=0):
    pad = "  " * depth
    cn = txt(rb, "ClassName") or "?"
    nm = txt(rb, "Name"); cid = txt(rb, "ID")
    print(f"{pad}[{rb.tag}] {nm}  class={cn.split('.')[-1] if cn else '?'} id={cid}")
    if rb.tag == "RigidBody":
        print(f"{pad}   mass={trip(rb,'Mass')} localPos={trip(rb,'LocalPosition')} rot={trip(rb,'Rotation')}")
        sz = rb.find("Size")
        if sz is not None:
            print(f"{pad}   size={trip(sz,'Length')} x {trip(sz,'Width')} x {trip(sz,'Height')}")
        # joint child
        for jt in rb.findall("Joint"):
            print(f"{pad}   JOINT {txt(jt,'Name')} id={txt(jt,'ID')} class={txt(jt,'ClassName')}")
            print(f"{pad}      localPos={trip(jt,'LocalPosition')} rot={trip(jt,'Rotation')} "
                  f"rotParent={trip(jt,'RotationRelativeParent')}")
            for c in jt:
                if c.tag in ("EnableLimits", "LowerLimit", "UpperLimit", "Limited", "MinLength", "MaxLength"):
                    print(f"{pad}      {c.tag}: {dict(c.attrib)}")
            for c in jt.findall("ConstraintRelaxation"):
                pass
    for ch in list(rb):
        if ch.tag in ("RigidBody", "Attachment"):
            dump_rb(ch, depth + 1)

dump_rb(body)
