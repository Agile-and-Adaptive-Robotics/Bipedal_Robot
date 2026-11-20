"""Convert the scaled subject01 OpenSim model using MyoConverter. Set GEOMETRY_ARG to an OpenSim geometry directory to override the bundled example geometry."""

import shutil
import sys
from pathlib import Path

import xml.etree.ElementTree as ET


def find_geometry():
    import os
    import myoconverter
    override = os.environ.get("GEOMETRY_ARG")
    if override:
        return Path(override)
    return Path(myoconverter.__file__).resolve().parents[1] / "models/osim/Gait2354Simbody/Geometry"


def find_repo_osim():
    return Path(__file__).resolve().with_name("subject01_simbody.osim")


def strip_missing_meshes(src, dst, geometry_folder):
    """Copy osim, dropping <Mesh> attached geometry whose file is missing
    (verbatim from the proven convert_to_mujoco.py)."""
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    root = ET.parse(src, parser=parser).getroot()
    dropped = 0
    for mesh in list(root.iter("Mesh")):
        mf = mesh.find("mesh_file")
        if mf is None or not (geometry_folder / (mf.text or "").strip()
                              ).exists():
            parent = next(p for p in root.iter() if mesh in list(p))
            parent.remove(mesh)
            dropped += 1
    ET.indent(ET.ElementTree(root), space="\t")
    ET.ElementTree(root).write(dst, encoding="UTF-8", xml_declaration=True)
    return dropped


if __name__ == "__main__":
    from myoconverter.O2MPipeline import O2MPipeline

    geometry = Path(sys.argv[1]) if len(sys.argv) > 1 else find_geometry()
    src = find_repo_osim()
    out = src.parent / "mjc" / "subject01_simbody"
    out.mkdir(parents=True, exist_ok=True)
    tmp_osim = out / "subject01_simbody.osim"
    n = strip_missing_meshes(src, tmp_osim, geometry)
    print(f"[subject01_simbody] stripped {n} missing-mesh decorations -> "
          f"{tmp_osim}", flush=True)
    print(f"[subject01_simbody] geometry: {geometry}", flush=True)
    O2MPipeline(
        str(tmp_osim), str(geometry), str(out),
        convert_steps=[1, 2, 3],
        muscle_list=None,
        osim_data_overwrite=True,
        conversion=True,
        validation=False,       # fast first pass (proven kwargs)
        speedy=True,
        generate_pdf=False,     # headless Windows
        add_ground_geom=True,
        treat_as_normal_path_point=False,
    )
    print(f"[subject01_simbody] DONE, see {out}", flush=True)
    print("NEXT: set AARL_MODEL to the new xml and run "
          "_validate_subject_model.py in Code/MuJoCo_SNS/spinal",
          flush=True)
