"""
Convert gait2392_robot.osim (Ben's mirrored robot model) and
gait2392_simbody.osim (pristine human model) to MuJoCo MJCF with MyoConverter.

Bone meshes (.vtp) are borrowed from myoconverter's Gait2354Simbody example,
which ships the identical OpenSim standard geometry set. Any mesh file that is
not available (e.g. treadmill.vtp) is stripped from a temporary conversion
copy of the model; MyoConverter re-adds a proper ground plane via
add_ground_geom=True.

Run inside the `myoconv` env:
    conda run -n myoconv python convert_to_mujoco.py

Conversion with validation/PDF off takes roughly 10-20 min per model on a
6-core laptop. Results land in ./mjc/<model_name>/.
"""

import shutil
from pathlib import Path

import xml.etree.ElementTree as ET

from myoconverter.O2MPipeline import O2MPipeline

HERE = Path(__file__).parent
GEOMETRY = Path(r"C:\Users\Ben\Documents\GitHub\myoconverter\models\osim\Gait2354Simbody\Geometry")
MODELS = ["gait2392_robot", "gait2392_simbody"]

kwargs = {
    "convert_steps": [1, 2, 3],
    "muscle_list": None,
    "osim_data_overwrite": True,
    "conversion": True,
    "validation": False,      # Vlt steps off: fast first pass; rerun later if needed
    "speedy": True,
    "generate_pdf": False,    # headless Windows: skip pyvista/GL report
    "add_ground_geom": True,
    "treat_as_normal_path_point": False,
}


def strip_missing_meshes(src, dst, geometry_folder):
    """Copy osim, dropping <Mesh> attached geometry whose file is missing."""
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    root = ET.parse(src, parser=parser).getroot()
    dropped = 0
    for mesh in list(root.iter("Mesh")):
        mf = mesh.find("mesh_file")
        if mf is None or not (geometry_folder / (mf.text or "").strip()).exists():
            parent = next(p for p in root.iter() if mesh in list(p))
            parent.remove(mesh)
            dropped += 1
    ET.indent(ET.ElementTree(root), space="\t")
    ET.ElementTree(root).write(dst, encoding="UTF-8", xml_declaration=True)
    return dropped


if __name__ == "__main__":
    for name in MODELS:
        out = HERE / "mjc" / name
        out.mkdir(parents=True, exist_ok=True)
        tmp_osim = out / f"{name}.osim"
        n = strip_missing_meshes(HERE / f"{name}.osim", tmp_osim, GEOMETRY)
        print(f"[{name}] stripped {n} missing-mesh decorations -> {tmp_osim}")
        O2MPipeline(str(tmp_osim), str(GEOMETRY), str(out), **kwargs)
        print(f"[{name}] DONE, see {out}")
