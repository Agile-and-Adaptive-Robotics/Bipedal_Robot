"""Build the two add_mag3_r P2 variant .osim copies from gait2392_robotbody.osim.

Ben's ruling (2026-09-26): the repo-vs-thumb-drive add_mag3_r P2 discrepancy is
settled by PERFORMANCE. This script writes two variant copies, touching ONLY
add_mag3_r-P2's pelvis-frame location:

  variant_repo  = (-0.059, -0.108, -0.03)   (as in repo gait2392_robotbody.osim;
                 identical numbers to add_mag3_r-P2_0)
  variant_thumb = (-0.1357, -0.0929, 0.0591) (full-precision original, master's
                 thumb drive D:/Bipedal humanoid/Gait2392_Robotbody/ -- not
                 mounted on easteregg2; coordinates from AGENTS.md 2026-09-22)

Originals are NOT modified.
"""

import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).parent
SRC = HERE.parent / "gait2392_robotbody.osim"
REPO_LOC = (-0.059, -0.108, -0.03)
THUMB_LOC = (-0.1357, -0.0929, 0.0591)


def load(path):
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    return ET.parse(path, parser=parser).getroot()


def find_p2(root):
    """Return the <location> element of add_mag3_r-P2 (the pelvis-frame one)."""
    for m in root.iter("Thelen2003Muscle"):
        if m.get("name") != "add_mag3_r":
            continue
        for pp in m.iter("PathPoint"):
            if pp.get("name") == "add_mag3_r-P2":
                frame = pp.find("socket_parent_frame").text.strip()
                loc = pp.find("location")
                assert "pelvis" in frame, f"P2 frame is {frame}, expected pelvis"
                return loc
    raise KeyError("add_mag3_r-P2 not found")


for tag, loc in (("repo", REPO_LOC), ("thumb", THUMB_LOC)):
    root = load(SRC)
    loc_el = find_p2(root)
    before = loc_el.text
    loc_el.text = f"{loc[0]!r} {loc[1]!r} {loc[2]!r}"
    out = HERE / f"variant_{tag}.osim"
    tree = ET.ElementTree(root)
    ET.indent(tree, space="\t")
    tree.write(out, encoding="UTF-8", xml_declaration=True)
    print(f"{out.name}: P2 {before!r} -> {loc_el.text!r}")

# Verify the two variants differ ONLY in that one element (and the repo variant
# is coordinate-identical to the source model).
r_src = load(SRC)
p_src = find_p2(r_src).text
r_repo = load(HERE / "variant_repo.osim")
r_thumb = load(HERE / "variant_thumb.osim")

def canon(root):
    """All path-point locations in document order, as canonical floats."""
    out = []
    for m in root.iter("Thelen2003Muscle"):
        for pp in m.iter():
            if pp.tag in ("PathPoint", "ConditionalPathPoint", "MovingPathPoint"):
                l = pp.find("location")
                v = tuple(float(x) for x in l.text.split()) if l is not None and l.text else None
                out.append((m.get("name"), pp.get("name"), v))
    return out

cs, crepo, cthumb = canon(r_src), canon(r_repo), canon(r_thumb)
assert cs == crepo, "repo variant unexpectedly differs from source"
diff = [(entry[0], entry[1], entry[2]) for entry, c in zip(cs, cthumb)
        if entry != c]
print(f"repo variant == source everywhere: OK ({len(cs)} path points)")
print(f"thumb variant differs in {len(diff)} point(s): {diff}")
assert len(diff) == 1 and diff[0][1] == "add_mag3_r-P2", "unexpected extra diffs"
print("VARIANT BUILD OK")
