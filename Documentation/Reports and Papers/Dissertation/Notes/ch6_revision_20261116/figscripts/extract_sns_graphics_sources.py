"""Read current SLX masks and six-edge reflex pathways without modifying models."""
from pathlib import Path
import hashlib
import json
import xml.etree.ElementTree as ET
import zipfile

root = next(p for p in Path(__file__).resolve().parents if p.name == "Bipedal_Robot")
src = root / "Code/Matlab/SNS_Simscape"
out = Path(__file__).parent

def system(path):
    with zipfile.ZipFile(path) as z:
        return ET.fromstring(z.read("simulink/systems/system_root.xml"))

icons = []
for b in system(src / "SNS_Library.slx").findall("Block"):
    mask = b.find("Mask")
    if mask is not None:
        icons.append(dict(name=b.attrib["Name"], display=mask.findtext("Display", ""),
                          description=mask.findtext("Description", "")))

models = {}
for name in ("KneeReflexDemo", "BeerCupReflexDemo"):
    tree = system(src / "demos" / (name + ".slx"))
    blocks = {b.attrib["SID"]: b for b in tree.findall("Block")}
    edges = []
    for line in tree.findall("Line"):
        source = line.findtext("P[@Name='Src']")
        if source:
            for dest in line.findall(".//P[@Name='Dst']"):
                edges.append((source.split("#")[0], dest.text.split("#")[0]))
    pathways = []
    for sid, b in blocks.items():
        if not b.attrib["Name"].startswith("syn_"):
            continue
        incoming = [blocks[a].attrib["Name"] for a, c in edges if c == sid]
        outgoing = [blocks[c].attrib["Name"] for a, c in edges if a == sid]
        p = {v.attrib["Name"]: v.text for v in b.findall("InstanceData/P")}
        pathways.append(dict(synapse=b.attrib["Name"], source=incoming, target=outgoing,
                             Esyn=p.get("Esyn"), gmax=p.get("gmax")))
    assert len(pathways) == 6, (name, pathways)
    assert sum(p["Esyn"] == "0" for p in pathways) == 2
    assert sum(p["Esyn"] == "-72" for p in pathways) == 4
    models[name] = pathways

files = [src / "SNS_Library.slx", src / "sns_icon_geom.m", src / "sns_icon_muscle.m",
         src / "demos/sns_build_demo.m", src / "demos/sns_build_beer_demo.m",
         src / "demos/KneeReflexDemo.slx", src / "demos/BeerCupReflexDemo.slx",
         src / "results/sns_beer_results.mat",
         src / "results/knee_transient_reflex_20261028/replay_summary.mat"]
hashes = {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
(out / "sns_graphics_sources.json").write_text(json.dumps(dict(icons=icons, pathways=models, sha256=hashes), indent=2), encoding="utf-8")
print(json.dumps(dict(blocks=len(icons), pathways=models), indent=2))
