"""smoke_test.py — verify headless rendering + helpers on this machine."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import render_common as rc  # noqa: E402
import bpy  # noqa: E402

scene, engine = rc.setup_scene(res=(480, 480))
print("ENGINE:", engine)
rc.add_sun()

root = rc.HERE
bone_mat = rc.make_material("bone", rc.BONE)
femur = rc.import_stl(rc.GEOM_DIR / "femur_r.stl", "femur", material=bone_mat)
print("FEMUR DIMS:", tuple(round(v, 3) for v in femur.dimensions))

cable = rc.poly_curve("cable", [(0, 0, 0), (0.2, 0.1, 0), (0.3, -0.05, 0)], 0.01,
                      material=rc.make_material("c", rc.CABLE_BLUE))
rc.keyframe_polyline(cable, [1, 20], [
    [(0, 0, 0), (0.2, 0.1, 0), (0.3, -0.05, 0)],
    [(0, 0, 0), (0.2, -0.1, 0), (0.3, 0.05, 0)],
])
lab = rc.text_label("Test label", "lbl", 0.05)
lab.location = (0.1, 0.2, 0.05)
cam = rc.look_at((0.1, -0.3, 0.6), (0.1, -0.2, 0.0), up=(0, 1, 0), lens=60)
rc.track_camera(lab, cam)
(root / "build").mkdir(exist_ok=True)
rc.render_frame(scene, 1, str(root / "build" / "smoke.png"))
rc.save_blend(root / "scenes" / "smoke.blend")
print("SMOKE_DONE engine=%s" % engine)
