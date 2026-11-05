"""build_visual3_frames.py — Gait2392 knee coordinate frames (right knee).

Kinematic reference visual: femur and tibia of the stock OpenSim Gait2392
model (gait2392_simbody.osim), muscles hidden, prescribed knee flexion.
Parent/child joint frame triads (femur_r_offset / tibia_r_offset, both at
zero offset from their bodies) move exactly per the CustomJoint:
    tibia origin in femur frame: p(q) = (tx(q), ty(q), 0)   [SimmSpline tables]
    tibia orientation: Rz(q)
Flexion is NEGATIVE (OpenSim convention, coordinate range -2.0944..0.1745 rad).
The same composition was verified against the converted MuJoCo model
(verify_gait2392_fk_20261104.py: 0.000 deg rotation error, <=0.87 mm position).

Run (background):
  blender --background --factory-startup --python build_visual3_frames.py -- \
      [--no-frames] [--stills-only]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import argparse
import math
import numpy as np
import bpy
from mathutils import Matrix, Vector
import render_common as rc

FPS = 15
DUR = 9.0                       # seconds (user asked 6-10 s slow loop)
N_FRAMES = int(DUR * FPS) + 1   # 136 -> loop seam frame N == frame 0
Q_MIN = -2.0944                 # rad, deep flexion (model range floor)
Q_MAX = 0.0                     # rad, straight (default pose)


def smoothstep(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3.0 - 2.0 * u)


def q_of_frame(k):
    """0 -> -2.0944 -> 0 ping-pong with eased turnarounds (flexion NEGATIVE)."""
    u = k / (N_FRAMES - 1)
    s = 1.0 - abs(2.0 * u - 1.0)          # 0->1->0 triangle
    return Q_MIN * smoothstep(s)


def interp_table(qgrid, values, q):
    return float(np.interp(q, qgrid, values))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-frames", action="store_true", help="build scene only")
    ap.add_argument("--stills-only", action="store_true")
    _argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    args, _ = ap.parse_known_args(_argv)

    scene, engine = rc.setup_scene(res=(720, 720), samples=32)
    rc.add_sun()
    print("ENGINE:", engine)

    gz = np.load(rc.HERE / "data" / "visual3_geom.npz")
    qgrid, txg, tyg = gz["q_grid"], gz["tx"], gz["ty"]

    bone_mat = rc.make_material("bone", (0.55, 0.54, 0.52))
    bone2_mat = rc.make_material("bone2", (0.47, 0.46, 0.44))
    femur = rc.import_stl(rc.GEOM_DIR / "femur_r.stl", "femur_r", material=bone_mat)
    tibia = rc.import_stl(rc.GEOM_DIR / "tibia_r.stl", "tibia_r", material=bone2_mat)

    frames = list(range(0, N_FRAMES))
    qs = [q_of_frame(k) for k in frames]
    locs = [(interp_table(qgrid, txg, q), interp_table(qgrid, tyg, q), 0.0) for q in qs]
    rots = [(0.0, 0.0, q) for q in qs]
    rc.keyframe_linear(tibia, "location", frames, locs)
    for i in range(3):
        rc.keyframe_linear(tibia, "rotation_euler", frames, [r[i] for r in rots], index=i)

    # ---- frame triads (parent static at hip; child rides the tibia) ------
    def triad(prefix, length=0.10, radius=0.007):
        parts = []
        for (vec, col, letter) in (((1, 0, 0), rc.AXIS_X, "x"),
                                   ((0, 1, 0), rc.AXIS_Y, "y"),
                                   ((0, 0, 1), rc.AXIS_Z, "z")):
            shaft, head = rc.arrow(f"{prefix}_{letter}", (0, 0, 0),
                                   tuple(np.array(vec) * length), radius,
                                   rc.make_material(f"{prefix}_{letter}_m", col))
            parts.append((shaft, head))
            tip = np.array(vec) * (length + 0.042)
            if letter == "z":
                # keep the letter clear of the arrowhead in side projection
                tip = tip + np.array((0.045, -0.045, 0.0))
            lab = rc.text_label(letter, f"{prefix}_{letter}_lab", 0.055, color=col)
            lab.location = tuple(tip)
            parts.append((lab, None))
        return parts

    femur_triad = triad("femframe")
    for ob, _ in femur_triad:
        pass  # world-placed, static (hip = femur frame origin)

    tibia_triad = triad("tibframe")
    for ob, _ in tibia_triad:
        ob.parent = tibia
        ob.matrix_parent_inverse.identity()

    # ---- labels (flat in the sagittal plane; readable from the +z side) ---
    cam = rc.look_at((-0.16, -0.42, 1.35), (-0.16, -0.40, 0.0), up=(0, 1, 0), lens=55)

    lab_femur = rc.text_label("Femur  (parent frame)", "lab_femur", 0.055,
                              align="CENTER")
    lab_femur.location = (-0.09, -0.155, 0.03)

    lab_tibia = rc.text_label("Tibia  (child frame)", "lab_tibia", 0.055,
                              align="CENTER")
    tib_lab_locs = [(tx - 0.10, ty + 0.075, 0.03) for (tx, ty, _) in locs]
    rc.keyframe_linear(lab_tibia, "location", frames, tib_lab_locs)

    lab_readout = rc.text_label("Knee angle:   0.0 deg", "lab_readout", 0.052,
                                align="CENTER")

    def update_readout(k):
        q = q_of_frame(k)
        deg = round(math.degrees(q))
        lab_readout.data.body = f"Knee angle: {deg:+d} deg" if deg != 0 else "Knee angle: 0 deg"

    # readout placement: fixed in CAMERA space -> child of the camera; the
    # text keeps its readable (+Z) side toward the camera view axis
    holder = bpy.data.objects.new("readout_holder", None)
    bpy.context.scene.collection.objects.link(holder)
    holder.parent = cam
    holder.matrix_parent_inverse.identity()
    holder.location = (0.0, 0.30, -1.0)
    lab_readout.parent = holder
    lab_readout.matrix_parent_inverse.identity()
    lab_readout.rotation_euler = (0, 0, 0)   # +Z_cam faces the viewer: readable

    rc.save_blend(rc.HERE / "scenes" / "visual3_knee_frames.blend")

    if args.stills_only:
        q = -1.6
        k = int(np.argmin([abs(q_of_frame(i) - q) for i in range(N_FRAMES)]))
        scene.render.resolution_x = 1080
        scene.render.resolution_y = 1080
        update_readout(k)
        scene.frame_set(k)
        rc.render_frame(scene, k, str(rc.HERE / "build" / "visual3" / "still_side.png"))
        # oblique companion: orbit the side camera 35 deg about the vertical
        target = Vector((-0.16, -0.40, 0.0))
        rel = Vector(cam.location) - target
        cam.location = target + Matrix.Rotation(math.radians(35), 3, "Y") @ rel
        cam.rotation_euler = rc.look_rotation(cam.location, target)
        cam.data.lens = 48
        rc.render_frame(scene, k, str(rc.HERE / "build" / "visual3" / "still_oblique.png"))
        print("STILLS_DONE")
        return

    if not args.no_frames:
        outdir = rc.HERE / "build" / "visual3" / "frames"
        outdir.mkdir(parents=True, exist_ok=True)
        for k in frames:
            update_readout(k)
            rc.render_frame(scene, k, str(outdir / f"f{k:04d}.png"))
            if k % 30 == 0:
                print(f"rendered {k}/{N_FRAMES}")
    print("VISUAL3_DONE")


if __name__ == "__main__":
    main()
