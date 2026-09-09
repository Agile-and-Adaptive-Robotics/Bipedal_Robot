"""
Add BPA tendon routes (site paths + zero-gain actuators) to a converted
MyoConverter MJCF model, or to any MuJoCo model.

The converted gait2392 model keeps its 92 Hill-type muscle actuators; this
helper appends NEW tendon routes for the BPAs that replace selected human
muscles (e.g. a 2x20 mm BPA set for the knee flexor). Each BPA gets:

    <site> nodes inside the named bodies defining the cable route
    <tendon><fixed>  over those sites
    <general> actuator with zero gain/bias  (force comes from BPAMuscleSystem)

Editing route points is then just editing the site pos="x y z" attributes in
plain XML - no OpenSim round-trip needed.

Usage
-----
    add_bpa(xml_in, xml_out, name="knee_flex_bpa_r",
            route=[("femur_r",  [0.006, -0.005, 0.055]),   # origin (hip frame side)
                   ("tibia_r",  [0.005, -0.033, 0.020])],  # insertion
            )
"""

import xml.etree.ElementTree as ET
from pathlib import Path


def _find_body(root, name):
    for b in root.iter("body"):
        if b.get("name") == name:
            return b
    raise KeyError(f"body '{name}' not found in MJCF")


def add_bpa(xml_in, xml_out, name, route, gear="0 0 0", ctrlrange="0 1"):
    """Append one BPA tendon + actuator. route = [(body, [x,y,z]), ...]."""
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    tree = ET.parse(xml_in, parser=parser)
    root = tree.getroot()

    site_names = []
    for i, (body_name, pos) in enumerate(route):
        body = _find_body(root, body_name)
        site = ET.SubElement(body, "site")
        site.set("name", f"{name}_s{i}")
        site.set("pos", " ".join(repr(float(v)) for v in pos))
        site.set("size", "0.005")
        site.set("rgba", "1 0.2 0.2 0.6")
        site_names.append(site.get("name"))

    tendon = ET.SubElement(root.find("tendon"), "spatial")
    tendon.set("name", f"{name}_tendon")
    for sn in site_names:
        s = ET.SubElement(tendon, "site")
        s.set("site", sn)

    act = ET.SubElement(root.find("actuator"), "general")
    act.set("name", name)
    act.set("tendon", f"{name}_tendon")
    act.set("gaintype", "fixed")
    act.set("gainprm", "0 0 0")
    act.set("biastype", "none")
    act.set("biasprm", "0 0 0")
    act.set("ctrlrange", ctrlrange)

    ET.indent(tree, space="  ")
    tree.write(xml_out, encoding="UTF-8", xml_declaration=True)
    return xml_out


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("xml_in")
    p.add_argument("xml_out")
    p.add_argument("name")
    p.add_argument("route", nargs="+",
                   help="body=x,y,z route points, e.g. femur_r=0.01,-0.02,0.05")
    a = p.parse_args()
    route = []
    for spec in a.route:
        body, xyz = spec.split("=", 1)
        route.append((body, [float(v) for v in xyz.split(",")]))
    print(add_bpa(a.xml_in, a.xml_out, a.name, route))
