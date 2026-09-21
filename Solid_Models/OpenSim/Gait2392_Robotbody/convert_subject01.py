"""Convert the SCALED subject01_simbody.osim to MuJoCo MJCF (MyoConverter).

Ben GO 2026-09-18: fixes the anthropometry mismatch found by
audit_ik_ground.py (Code/MuJoCo_SNS/spinal) — the IK/GRF reference data
is from the marker-scaled subject01 (standing pelvis ~1.02 m) while our
gait2392_simbody_cvt3.xml is the DEFAULT unscaled model (0.95 m), which
leaves the replayed stance foot ~8.7 cm above the floor.

RUN THIS ON EASTEREGG2 (the machine with the patched MyoConverter clone
and the opensim env):
    cd D:\GitHub\myoconverter
    copy /Y D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\convert_subject01.py .
    set CONDA_PREFIX=D:\Anaconda\envs\myo
    D:\Anaconda\envs\myo\python.exe convert_subject01.py
(the `if __name__ == "__main__"` guard is REQUIRED - MyoConverter step 3
multiprocessing-spawns children; a guard-less child re-runs the whole
pipeline and deadlocks on the log file.)

The script auto-detects clone layout for the Geometry folder:
  - easteregg2:  D:\GitHub\myoconverter\models\osim\Gait2354Simbody\Geometry
  - laptop:      C:\Users\Ben\Documents\GitHub\myoconverter\models\...
Bone meshes are the same standard Gait2354/Gait2392 set (scaling does
not change mesh files); missing decorations (treadmill.vtp) are stripped
exactly as convert_to_mujoco.py does.

Output: mjc/subject01_simbody/ next to this script's repo folder.
Expected wall time: 10-20 min (validation/PDF off, like the proven run).

AFTER IT FINISHES (on any machine, myo env):
    cd Code\MuJoCo_SNS\spinal
    set AARL_MODEL=<repo>\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\subject01_simbody\<...>.xml
    python _validate_subject_model.py
which checks standing height ~1.02 m, IK-replay stance foot ON the
floor, 92 actuators, and finite dynamics - the gates that motivated the
conversion.  Do NOT wire the new model into tuning until the audit
report says PASS.
"""
import shutil
import sys
from pathlib import Path

import xml.etree.ElementTree as ET


def find_geometry():
    here = Path(__file__).resolve().parent
    cands = [
        here / "models" / "osim" / "Gait2354Simbody" / "Geometry",
        Path(r"D:\GitHub\myoconverter\models\osim\Gait2354Simbody\Geometry"),
        Path(r"C:\Users\Ben\Documents\GitHub\myoconverter\models\osim"
             r"\Gait2354Simbody\Geometry"),
    ]
    for c in cands:
        if c.is_dir():
            return c
    raise FileNotFoundError(
        "Gait2354Simbody Geometry folder not found; pass it via "
        "GEOMETRY_ARG env var")


def find_repo_osim():
    here = Path(__file__).resolve().parent
    cands = [
        here / "subject01_simbody.osim",
        Path(r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim"
             r"\Gait2392_Robotbody\subject01_simbody.osim"),
        Path(r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim"
             r"\Gait2392_Robotbody\subject01_simbody.osim"),
    ]
    for c in cands:
        if c.is_file():
            return c
    raise FileNotFoundError("subject01_simbody.osim not found")


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
