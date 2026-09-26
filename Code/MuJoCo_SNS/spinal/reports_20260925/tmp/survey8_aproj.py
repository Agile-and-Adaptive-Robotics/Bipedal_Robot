# -*- coding: utf-8 -*-
"""Dump all joints (limits), springs, and one full muscle."""
import xml.etree.ElementTree as ET

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()

def txt(el, tag):
    c = el.find(tag)
    return (c.text or "").strip() if c is not None else None

def sub(el, tag):
    return el.find(tag) if el is not None else None

def triplet(el, tag, attr="Actual"):
    c = el.find(tag) if el is not None else None
    return c.get(attr) if c is not None else None

org = root.find(".//Organism")

print("=== ALL JOINTS ===")
for jt in org.iter("Joint"):
    nm = txt(jt, "Name")
    lo = sub(jt, "LowerLimit"); hi = sub(jt, "UpperLimit")
    print(f"{nm}: EnableLimits={txt(jt,'EnableLimits')} "
          f"lower={triplet(lo,'LimitPos')} upper={triplet(hi,'LimitPos')} "
          f"restLo={triplet(lo,'Restitution')} "
          f"MaxForce={triplet(jt,'MaxForce')} motor={txt(jt,'EnableMotor')} type={txt(jt,'MotorType')}")

print("\n=== SPRINGS ===")
# Spring is a RigidBody Type=Spring
for rb in org.iter("RigidBody"):
    if txt(rb, "Type") == "Spring":
        print("Name:", txt(rb, "Name"), "ID:", txt(rb, "ID"))
        print("  children:", [c.tag for c in rb])
        for c in rb:
            if c.tag not in ("Name", "ID", "Description", "Transparencies", "IsVisible", "Ambient",
                             "Diffuse", "Specular", "Shininess", "Texture", "LocalPosition",
                             "Rotation", "LocalMatrix", "DraggerSize", "Type", "PartType"):
                print("  ", c.tag, (dict(c.attrib) if c.attrib else (c.text or "")[:80]))
                for cc in c:
                    print("      ", cc.tag, dict(cc.attrib) if cc.attrib else (cc.text or "")[:60])

print("\n=== FULL MUSCLE hip_L_flx ===")
for rb in org.iter("RigidBody"):
    if txt(rb, "Type") == "LinearHillMuscle" and txt(rb, "Name") == "hip_L_flx":
        for c in rb:
            attrs = dict(c.attrib) if c.attrib else ""
            text = (c.text or "").strip()
            print(f"<{c.tag}> {attrs} {text[:90]}")
            for cc in c:
                print(f"   <{cc.tag}> {dict(cc.attrib) if cc.attrib else (cc.text or '').strip()[:60]}")
        break

print("\n=== MUSCLE ApplyTension check ===")
for rb in org.iter("RigidBody"):
    if txt(rb, "Type") == "LinearHillMuscle":
        print(f"  {txt(rb,'Name')}: ApplyTension={txt(rb,'ApplyTension')}")
