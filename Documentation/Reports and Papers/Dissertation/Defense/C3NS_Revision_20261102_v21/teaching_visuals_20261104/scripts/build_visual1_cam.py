"""build_visual1_cam.py — cam-assisted extensor teaching visual.

Every rendered quantity comes from the ADOPTED cam design of record
(Vas_Pam_20mm_Result_cam_20261014_0416.mat -> data/visual1_geom.npz):
  route p1..p9 per pose (rows 1-5 femur frame, 6-9 t1 -> femur via the
  evaluator's T_ICR_t1 / T_Pam transforms), the migrating knee pivot
  (ICR, T_Pam translation), the tibial clip-part cam surface arc
  (circumcircle of the full-flexion contacts p6/p7/p8; the 96.9 mm value
  of record is the loss-model drum radius p7x-p9x, noted in the README),
  the femoral cylinder pulley (geo.femurCylCenter/radius), and the
  moment arm of record (perpendicular from the pivot to the p6 line of
  action; validated against predBest.momentArm, mean diff 0.18 mm).
Bones: OpenSim femur_r/tibia_r STLs registered exactly as in the
dissertation figures (appendix_c_bone_background.m offsets).

Run:
  blender --background --factory-startup --python build_visual1_cam.py -- [--stills-only]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import argparse
import math
import numpy as np
import bpy
from mathutils import Vector, Matrix
import render_common as rc

FPS = 15
INTRO = 18                      # frames 0..17   oblique -> side + bone fade
SWEEP_A = 42                    # frames 18..59  flexion -> extension
SWEEP_B = 42                    # frames 60..101 extension -> flexion
OUTRO = 18                      # frames 102..119 side -> oblique + unfade
N_FRAMES = INTRO + SWEEP_A + SWEEP_B + OUTRO   # 120 = 8 s at 15 fps
BONE_ALPHA_ON = 1.0
BONE_ALPHA_GHOST = 0.22


def smoothstep(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3.0 - 2.0 * u)


def build_schedule(phiD):
    """Frame list -> phi (deg), with held flexion during intro/outro."""
    phi_min, phi_max = float(phiD.min()), float(phiD.max())   # -120, +10
    phis = np.zeros(N_FRAMES)
    phis[:INTRO] = phi_min
    u = (np.arange(SWEEP_A) + 1) / SWEEP_A
    phis[INTRO:INTRO + SWEEP_A] = phi_min + (phi_max - phi_min) * u
    u = np.arange(SWEEP_B) / (SWEEP_B - 1)
    phis[INTRO + SWEEP_A:INTRO + SWEEP_A + SWEEP_B] = phi_max + (phi_min - phi_max) * u
    phis[INTRO + SWEEP_A + SWEEP_B:] = phi_min
    # ease the sweep ends
    phis[INTRO:INTRO + SWEEP_A] = (
        phi_min + (phi_max - phi_min) * smoothstep((np.arange(SWEEP_A) + 1) / SWEEP_A))
    seg = smoothstep(np.arange(SWEEP_B) / (SWEEP_B - 1))
    phis[INTRO + SWEEP_A:INTRO + SWEEP_A + SWEEP_B] = phi_max + (phi_min - phi_max) * seg
    return phis


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stills-only", action="store_true")
    _argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    args, _ = ap.parse_known_args(_argv)

    scene, engine = rc.setup_scene(res=(720, 720), samples=32)
    rc.add_sun()
    print("ENGINE:", engine)

    g = np.load(rc.HERE / "data" / "visual1_geom.npz")
    phiD = g["phiD"]
    route = g["route_f"]            # (9, 3, N) femur frame
    active = g["active"]
    pivot = g["pivot"]              # (N, 3)
    tib_rot = g["tibia_rot"]        # (N,) rad about z
    t1toICR = g["t1toICR"]          # (N, 3) t1-origin offset from the pivot
    loa_p, loa_d, foot, mA = g["loa_p"], g["loa_d"], g["foot"], g["moment_arm"]
    dc = g["drum_center_t1"]
    seeds = g["seeds_t1"]
    p_end = g["p_end_t1"]
    pul_c, pul_r = g["femur_pulley_center"], g["femur_pulley_radius"]
    f_off, t_off = g["femur_reg_off"], g["tibia_reg_off"]
    p1 = route[0]                   # (3, N) femur-fixed origin

    frames = list(range(0, N_FRAMES))
    phis = build_schedule(phiD)

    def phi_to_idx(phi):
        return int(np.argmin(np.abs(phiD - phi)))

    idxs = [phi_to_idx(p) for p in phis]

    # ---- materials ---------------------------------------------------------
    femur_mat = rc.make_material("femur_m", rc.BONE)
    tibia_mat = rc.make_material("tibia_m", (0.47, 0.46, 0.44))
    pulley_mat = rc.make_material("pulley_m", rc.PULLEY_GRAY)
    cable_mat = rc.make_material("cable_m", rc.CABLE_BLUE)
    arm_mat = rc.make_material("arm_m", rc.ARM_GOLD)

    # ---- bones (registered exactly like the dissertation figures) ---------
    femur = rc.import_stl(rc.GEOM_DIR / "femur_r.stl", "femur", offset=f_off,
                          material=femur_mat)
    tibia = rc.import_stl(rc.GEOM_DIR / "tibia_r.stl", "tibia", offset=t_off,
                          material=tibia_mat)
    # tibia object: local = t1-frame coords; world = Rz(phi)*(p - t1toICR) + pivot
    locs = []
    for i in idxs:
        c, s = math.cos(float(tib_rot[i])), math.sin(float(tib_rot[i]))
        Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        locs.append(tuple(pivot[i] - Rz @ t1toICR[i]))
    rots = [(0.0, 0.0, float(tib_rot[i])) for i in idxs]
    rc.keyframe_linear(tibia, "location", frames, locs)
    for i in range(3):
        rc.keyframe_linear(tibia, "rotation_euler", frames,
                           [r[i] for r in rots], index=i)

    # ---- shaped cam: the actual TibiaClip_002 hardware (tibia-fixed) -------
    # Registered to the evaluator t1 frame: wrap-circle center match (radius
    # 59.83 mm verified on the real mesh) + the p8 "clip-edge" endpoint rule
    # (see prep_geometry.py / README). The two halves straddle the sagittal
    # plane with the cable slot at z=0.
    clip_reg = np.load(rc.HERE / "data" / "clip_registration.npz")
    C_P = clip_reg["C_P"]                      # mm, part frame
    delta = math.radians(float(clip_reg["delta_deg"]))
    C_T1 = clip_reg["C_T1"]                    # mm, t1 frame
    z_mid = float(clip_reg["z_mid"])
    c, s = math.cos(delta), math.sin(delta)
    R = np.array([[c, -s], [s, c]])

    def clip_transform(verts_mm):
        xy = verts_mm[:, :2] - C_P
        xy = xy @ R.T + C_T1
        z = verts_mm[:, 2] - z_mid
        return np.column_stack([xy, z]) / 1000.0

    cam_mat = rc.make_material("cam_hw", (0.50, 0.53, 0.57))   # machined steel
    for side in ("lateral", "medial"):
        stl_path = str(clip_reg[f"stl_{side}"])
        ob = rc.import_stl_transformed(stl_path, f"clip_{side}",
                                       clip_transform, cam_mat)
        ob.parent = tibia

    # contact markers (tibia-fixed)
    for i in range(3):
        bpy.ops.mesh.primitive_uv_sphere_add(
            radius=0.006, location=(seeds[i, 0], seeds[i, 1], 0.0))
        m = bpy.context.active_object
        m.name = f"contact{i}"
        rc.set_material(m, rc.make_material(f"contact{i}_m", (0.05, 0.05, 0.06)))
        m.parent = tibia

    # insertion marker (tibia-fixed)
    bpy.ops.mesh.primitive_uv_sphere_add(
        radius=0.008, location=(p_end[0], p_end[1], 0.0))
    ins = bpy.context.active_object
    ins.name = "insertion"
    rc.set_material(ins, rc.make_material("ins_m", rc.CABLE_BLUE))
    ins.parent = tibia

    # femoral pulley (parent-frame fixed) + origin marker
    rc.cylinder("femur_pulley", pul_r, 0.02, (pul_c[0], pul_c[1], 0.0), pulley_mat)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.008,
                                         location=(p1[:, 0][0], p1[:, 1][0], 0.0))
    org = bpy.context.active_object
    org.name = "origin"
    rc.set_material(org, rc.make_material("org_m", rc.ARM_GOLD))

    # ---- force path (the route of record, one tube; the dissertation figure
    # draws the whole optimized route as the "BPA path") ---------------------
    route_pts = [ [tuple(route[j, :, i]) for j in range(9)] for i in idxs ]
    cable = rc.poly_curve("cable", route_pts[0], radius=0.005, material=cable_mat)
    rc.keyframe_polyline(cable, frames, route_pts)

    # ---- moment arm + pivot marker -----------------------------------------
    pivot_marker = bpy.data.meshes.new("pivot_marker")
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.009, location=tuple(pivot[idxs[0]]))
    pm = bpy.context.active_object
    pm.name = "pivot_marker"
    rc.set_material(pm, rc.make_material("pm_m", rc.ARM_GOLD))
    arm_locs = [tuple(pivot[i]) for i in idxs]
    rc.keyframe_linear(pm, "location", frames, arm_locs)

    arm_pts = [[tuple(pivot[i]), tuple(foot[i])] for i in idxs]
    arm = rc.poly_curve("moment_arm", arm_pts[0], radius=0.005, material=arm_mat)
    rc.keyframe_polyline(arm, frames, arm_pts)

    # ---- force arrow along the line of action ------------------------------
    arrow_pts = []
    head_locs, head_rots = [], []
    for i in idxs:
        p6 = Vector(loa_p[i])
        d = Vector(loa_d[i])            # unit, p6 -> p5 (pull direction)
        head = p6 + d * 0.09
        arrow_pts.append([tuple(p6), tuple(head - d * 0.025)])
        head_locs.append(tuple(head))
        head_rots.append(d.to_track_quat("Z", "Y").to_euler())
    f_shaft = rc.poly_curve("force_shaft", arrow_pts[0], radius=0.006, material=arm_mat)
    rc.keyframe_polyline(f_shaft, frames, arrow_pts)
    cone_mesh = bpy.data.meshes.new("force_head")
    bpy.ops.mesh.primitive_cone_add(radius1=0.016, depth=0.032, vertices=24,
                                    location=(0, 0, 0))
    cone_ob = bpy.context.active_object
    cone_mesh = cone_ob.data
    cone_mesh.name = "force_headmesh"
    bpy.data.objects.remove(cone_ob)
    f_head = bpy.data.objects.new("force_head", cone_mesh)
    bpy.context.scene.collection.objects.link(f_head)
    rc.set_material(f_head, arm_mat)
    rc.keyframe_linear(f_head, "location", frames, head_locs)
    for i in range(3):
        rc.keyframe_linear(f_head, "rotation_euler", frames,
                           [r[i] for r in head_rots], index=i)

    # ---- labels -------------------------------------------------------------
    texts = []
    leaders = []

    def add_text(body, name, size=0.028, color=rc.INK, align="CENTER"):
        ob = rc.text_label(body, name, size, color=color, align=align)
        ob.location = (ob.location.x, ob.location.y, 0.012)
        texts.append(ob)
        return ob

    def add_leader(name, a, b, radius=0.0015):
        ob = rc.poly_curve(name, [a, b], radius=radius,
                           material=rc.make_material(name + "_m", (0.35, 0.35, 0.38)))
        leaders.append(ob)
        return ob

    lab_pivot = add_text("Knee pivot", "L_pivot")
    lab_cam = add_text("Shaped cam", "L_cam")
    lab_bpa = add_text("BPA / cable", "L_bpa", color=rc.CABLE_BLUE)
    lab_arm = add_text("Moment arm", "L_arm", color=rc.ARM_GOLD)
    lab_force = add_text("BPA force", "L_force", color=rc.ARM_GOLD)

    # label/leader keyframes
    piv_lab_locs, piv_leaders = [], []
    cam_lab_locs, cam_leaders = [], []
    arm_lab_locs, arm_leaders = [], []
    force_lab_locs, force_leaders = [], []
    for i in idxs:
        pv = Vector(pivot[i])
        piv_lab_locs.append(tuple(pv + Vector((-0.085, 0.052, 0.012))))
        piv_leaders.append([tuple(pv + Vector((-0.012, 0.010, 0.004))),
                            tuple(pv + Vector((-0.072, 0.044, 0.004)))])
        dcv = Vector((dc[0], dc[1], 0.0))
        r_cam = float(g["drum_radius"])
        Rz = Matrix.Rotation(tib_rot[i], 3, "Z")
        dcv_w = Rz @ (dcv - Vector(t1toICR[i]))   # t1 -> femur via the tibia object
        cam_lab_locs.append(tuple(Vector(pivot[i]) + dcv_w + Vector((0.075, 0.012, 0.012))))
        cam_leaders.append([tuple(Vector(pivot[i]) + dcv_w + Vector((r_cam * 0.45, -r_cam * 0.45, 0.004))),
                            tuple(Vector(pivot[i]) + dcv_w + Vector((0.060, -0.004, 0.004)))])
        ft = Vector(foot[i])
        arm_lab_locs.append(tuple(ft + Vector((-0.020, -0.090, 0.012))))
        arm_leaders.append([tuple(ft + Vector((0.0, -0.006, 0.004))),
                            tuple(ft + Vector((-0.016, -0.078, 0.004)))])
        p6 = Vector(loa_p[i])
        d = Vector(loa_d[i])
        hv = p6 + d * 0.09
        force_lab_locs.append(tuple(hv + Vector((0.095, 0.048, 0.012))))
        force_leaders.append([tuple(hv + Vector((0.012, 0.004, 0.004))),
                              tuple(hv + Vector((0.082, 0.040, 0.004)))])
    rc.keyframe_linear(lab_pivot, "location", frames, piv_lab_locs)
    rc.keyframe_linear(lab_cam, "location", frames, cam_lab_locs)
    rc.keyframe_linear(lab_arm, "location", frames, arm_lab_locs)
    rc.keyframe_linear(lab_force, "location", frames, force_lab_locs)
    lab_bpa.location = (float(p1[0, 0]) - 0.055, float(p1[1, 0]) - 0.045, 0.012)
    lead_bpa = add_leader("lead_bpa",
                          (float(p1[0, 0]) - 0.010, float(p1[1, 0]) - 0.006, 0.004),
                          (float(p1[0, 0]) - 0.048, float(p1[1, 0]) - 0.038, 0.004))
    # "Moment arm" / "BPA force" texts fade in once the arm is readable
    # (deep flexion hides it behind the cam) and fade out before the outro
    fade_pairs = []
    for ob in (lab_arm, lab_force):
        fade_pairs.append(ob.data.materials[0])
    for ob in leaders:
        if ob.name in ("lead_arm", "lead_force"):
            fade_pairs.append(ob.data.materials[0])
    fade_frames, fade_vals = [], []
    for k in frames:
        if k < 24:
            v = 0.0
        elif k < 30:
            v = (k - 24) / 6.0
        elif k <= 96:
            v = 1.0
        elif k < 102:
            v = 1.0 - (k - 96) / 6.0
        else:
            v = 0.0
        fade_frames.append(k)
        fade_vals.append(max(0.0, min(1.0, v)))
    for mat in fade_pairs:
        rc.keyframe_mat_alpha(mat, fade_frames, fade_vals)
    for ob, pts_list in (
            (add_leader("lead_pivot", piv_leaders[0][0], piv_leaders[0][1]), piv_leaders),
            (add_leader("lead_cam", cam_leaders[0][0], cam_leaders[0][1]), cam_leaders),
            (add_leader("lead_arm", arm_leaders[0][0], arm_leaders[0][1]), arm_leaders),
            (add_leader("lead_force", force_leaders[0][0], force_leaders[0][1]), force_leaders)):
        rc.keyframe_polyline(ob, frames, pts_list)

    # ---- camera: oblique intro -> side view teaching sweep ------------------
    # frame the whole registered assembly across every pose (the shank tip
    # reaches x ~ -0.34 / y ~ -0.79 at deep flexion and extension)
    ylo = -0.84
    yhi = float(route[0, 1, :].max() + 0.06)
    xlo = -0.40
    xhi = float(max(route[0, 0, :].max(), 0.20))
    target = Vector((0.5 * (xlo + xhi), 0.5 * (ylo + yhi), 0.0))
    half_h = 0.5 * (yhi - ylo) * 1.12
    dist = half_h / (18.0 / 50.0)                     # lens 50 mm, 36 mm sensor
    side_loc = Vector((target.x, target.y, dist))
    cam = rc.look_at(tuple(side_loc), tuple(target), lens=50)

    cam_locs, cam_rots = [], []
    for k in frames:
        # sample the intro/outro fades mid-frame so frame 0 and the loop's
        # last frame sit at the same point of the swing (seamless GIF)
        if k < INTRO:
            u = smoothstep((k + 0.5) / INTRO)
        elif k >= INTRO + SWEEP_A + SWEEP_B:
            u = 1.0 - smoothstep((k - (INTRO + SWEEP_A + SWEEP_B) + 0.5) / OUTRO)
        else:
            u = 1.0
        ang = math.radians(42 * (1.0 - u))
        rel = side_loc - target
        pos = target + Matrix.Rotation(-ang, 3, "Y") @ rel
        cam_locs.append(tuple(pos))
        cam_rots.append(rc.look_rotation(pos, target))
    rc.keyframe_linear(cam, "location", frames, cam_locs)
    for i in range(3):
        rc.keyframe_linear(cam, "rotation_euler", frames,
                           [r[i] for r in cam_rots], index=i)

    # readout (camera child)
    readout = rc.text_label("Knee angle:  -120 deg", "L_readout", 0.038, align="CENTER")
    holder = bpy.data.objects.new("readout_holder", None)
    bpy.context.scene.collection.objects.link(holder)
    holder.parent = cam
    holder.matrix_parent_inverse.identity()
    holder.location = (0.0, 0.30, -1.0)
    readout.parent = holder
    readout.matrix_parent_inverse.identity()

    def update_readout(k):
        readout.data.body = f"Knee angle: {phis[k]:+.0f} deg"

    # ---- bone fade -----------------------------------------------------------
    for mat in (femur_mat, tibia_mat):
        vals = []
        for k in frames:
            if k < INTRO:
                u = smoothstep((k + 0.5) / INTRO)
            elif k >= INTRO + SWEEP_A + SWEEP_B:
                u = 1.0 - smoothstep((k - (INTRO + SWEEP_A + SWEEP_B) + 0.5) / OUTRO)
            else:
                u = 1.0
            vals.append(BONE_ALPHA_ON + (BONE_ALPHA_GHOST - BONE_ALPHA_ON) * u)
        rc.keyframe_mat_alpha(mat, frames, vals)

    rc.save_blend(rc.HERE / "scenes" / "visual1_cam_extensor.blend")

    if args.stills_only:
        # labeled still at the first up-sweep frame at/after -30 deg
        sweep = phis[INTRO:INTRO + SWEEP_A]
        k = INTRO + int(np.argmax(sweep >= -30.0))
        scene.render.resolution_x = 1080
        scene.render.resolution_y = 1080
        update_readout(k)
        scene.frame_set(k)
        rc.render_frame(scene, k, str(rc.HERE / "build" / "visual1" / "still_labeled.png"))
        for o in texts:
            o.hide_render = True
        for o in leaders:
            o.hide_render = True
        rc.render_frame(scene, k, str(rc.HERE / "build" / "visual1" / "still_unlabeled.png"))
        print("STILLS_DONE")
        return

    outdir = rc.HERE / "build" / "visual1" / "frames"
    outdir.mkdir(parents=True, exist_ok=True)
    for k in frames:
        update_readout(k)
        rc.render_frame(scene, k, str(outdir / f"f{k:04d}.png"))
        if k % 30 == 0:
            print(f"rendered {k}/{N_FRAMES}")
    print("VISUAL1_DONE")


if __name__ == "__main__":
    main()
