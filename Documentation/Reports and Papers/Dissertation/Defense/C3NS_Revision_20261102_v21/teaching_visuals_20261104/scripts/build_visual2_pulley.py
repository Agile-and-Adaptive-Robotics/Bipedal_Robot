"""build_visual2_pulley.py — reverse tackle with two parallel BPAs (3D).

Arrangement of record: Figures/30-torque/reverse_pulley_biomimetic.pdf and
reverse_pulley_principle.pdf, generators
  Components/reverse_pulley_principle/draw_biomimetic_reverse_pulley_20261029.py
  Components/reverse_pulley_principle/draw_reverse_pulley_publication.py
All 2-D schematic coordinates below are taken verbatim from those scripts
(femur/tibia polygons, BPA spans, sheave centers, cable tangents). The
animation shows the IDEAL kinematics: each BPA contracts by delta, the moving
sheave translates delta, and the output tendon past the fixed redirect feeds
out by 2*delta (two line parts between moving sheave and fixed mount).
The ideal force relation (output = half the combined BPA force) is labeled;
measured efficiency and dynamic BPA validation are deliberately NOT shown.

Run:
  blender --background --factory-startup --python build_visual2_pulley.py -- [--stills-only]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import argparse
import math
import numpy as np
import bpy
import render_common as rc
import pulley_kinematics as pk

FPS = 15
DUR = 7.0
N_FRAMES = int(DUR * FPS) + 1
Z_TH = 0.12          # schematic plate thickness (display only)
DELTA_MAX = 0.42     # BPA contraction, schematic units
LABEL_Z = 0.10       # text plane in front of the hardware

# figure geometry (draw_biomiketic_reverse_pulley_20261029.py coordinates)
FEMUR_POLY = [[.6, 3.8], [.6, 4.5], [5.9, 4.6], [6.7, 4.4], [7.1, 3.8], [6.7, 3.1], [6.0, 3.5]]
TIBIA_POLY = [[6.9, 3.2], [7.45, 3.4], [8.75, .65], [8.15, .4], [6.7, 2.9]]
KNEE_DOTS = [(6.25, 3.85), (6.65, 4.05), (6.75, 3.2), (7.05, 3.35)]
KNEE_LINKS = [((6.25, 3.85), (7.05, 3.35)), ((6.65, 4.05), (6.75, 3.2))]
MOUNT_X = 0.85                     # origin mount bar
BPA_Y = (5.55, 6.45)               # the two parallel BPAs
BPA_X0, BPA_X1 = 1.15, 3.8         # BPA body span (at delta = 0)
STUB_END = 4.08                    # fixed-length stub after the BPA body
SHEAVE0 = (4.8, 6.0)               # moving sheave center at delta = 0
SHEAVE_R = 0.45
ANCHOR = (7.35, 6.45)              # cable anchor on the fixed tackle mount
REDIRECT = (7.35, 5.1)             # fixed redirect sheave center
TENDON_X = 7.8                     # output tendon line (tangent of redirect)
INSERTION0 = (TENDON_X, 2.35)      # insertion diamond at delta = 0


def smoothstep(u):
    u = np.clip(u, 0.0, 1.0)
    return u * u * (3.0 - 2.0 * u)


def delta_of_frame(k):
    u = k / (N_FRAMES - 1)
    return DELTA_MAX * smoothstep(1.0 - abs(2.0 * u - 1.0))


def plate(name, poly, material, z=Z_TH):
    verts3 = [(p[0], p[1], 0.0) for p in poly]
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts3, [], [list(range(len(poly)))])
    mesh.update()
    ob = bpy.data.objects.new(name, mesh)
    bpy.context.scene.collection.objects.link(ob)
    sol = ob.modifiers.new("Solidify", "SOLIDIFY")
    sol.thickness = z
    sol.offset = 0.0
    rc.set_material(ob, material)
    return ob


def cable_points(delta):
    knee = knee_pose(delta)
    sx = SHEAVE0[0] - delta                       # moving sheave center x
    pts = [(ANCHOR[0], ANCHOR[1], 0)]
    pts.append((sx, SHEAVE0[1] + SHEAVE_R, 0))    # top tangent of moving sheave
    for a in np.linspace(90, 270, 9)[1:-1]:       # half wrap (deg)
        pts.append((sx + SHEAVE_R * math.cos(math.radians(a)),
                    SHEAVE0[1] + SHEAVE_R * math.sin(math.radians(a)), 0))
    pts.append((sx, SHEAVE0[1] - SHEAVE_R, 0))    # bottom tangent
    pts.append((REDIRECT[0], REDIRECT[1] + SHEAVE_R, 0))   # top of redirect
    for a in np.linspace(90, math.degrees(knee['angle']), 6)[1:-1]:
        pts.append((REDIRECT[0] + SHEAVE_R * math.cos(math.radians(a)),
                    REDIRECT[1] + SHEAVE_R * math.sin(math.radians(a)), 0))
    pts.append((*knee['tangent'], 0))
    pts.append((*knee['P'], 0))
    return pts


def knee_pose(delta):
    lo, hi = 0.0, pk.end_angle(DELTA_MAX)
    for _ in range(40):
        mid=(lo+hi)/2
        if pk.contraction(mid)<delta:lo=mid
        else:hi=mid
    return pk.pose((lo+hi)/2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stills-only", action="store_true")
    _argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
    args, _ = ap.parse_known_args(_argv)

    scene, engine = rc.setup_scene(res=(1080, 720), samples=32)
    rc.add_sun(energy=0.4, direction=(-0.3, -0.6, 1.0))
    print("ENGINE:", engine)

    bone_m = rc.make_material("femur_m", (0.78, 0.75, 0.70))
    bone2_m = rc.make_material("tibia_m", (0.72, 0.69, 0.63))
    red_m = rc.make_material("bpa_m", rc.BPA_RED)
    grey_m = rc.make_material("grey_m", rc.PULLEY_GRAY)
    grey_d = rc.make_material("grey_d", (0.42, 0.43, 0.45))
    cable_m = rc.make_material("cable_m", rc.CABLE_BLUE)

    plate("femur", FEMUR_POLY, bone_m)
    plate("tibia", TIBIA_POLY, bone2_m)
    # biomimetic crossed knee linkage (static in this schematic)
    for i, ((x1, y1), (x2, y2)) in enumerate(KNEE_LINKS):
        seg = rc.poly_curve(f"kneelink{i}", [(x1, y1, 0), (x2, y2, 0)], 0.035,
                            material=grey_d)
        seg.location = (0, 0, Z_TH / 2)
    for i, (x, y) in enumerate(KNEE_DOTS):
        bpy.ops.mesh.primitive_uv_sphere_add(radius=0.06, location=(x, y, Z_TH / 2))
        d = bpy.context.active_object
        d.name = f"kneedot{i}"
        rc.set_material(d, rc.make_material(f"kneedot{i}_m", (0.95, 0.95, 0.95)))

    # origin mount bar
    mount = rc.box("origin_mount", (0.10, 2.4, Z_TH), (MOUNT_X, 5.7, 0), grey_d)
    for yy in BPA_Y:
        rc.box(f"mountstub{yy}", (0.30, 0.06, Z_TH), (MOUNT_X + 0.15, yy, 0), grey_d)

    # two parallel BPAs + their distal stubs + crosshead + moving sheave
    bpa_len0 = BPA_X1 - BPA_X0
    bpas = []
    for yy in BPA_Y:
        ob = rc.box(f"bpa{yy}", (bpa_len0, 0.46, Z_TH), (0, 0, 0), red_m)
        bpas.append(ob)
        rc.box(f"stub{yy}", (0.28, 0.06, Z_TH), (0, yy, 0), grey_d)
    cross = rc.box("crosshead", (0.06, SHEAVE0[1] - BPA_Y[0] + SHEAVE_R, Z_TH),
                   (0, 0, 0), grey_d)
    cross_arm = rc.box("crossarm", (SHEAVE0[0] - STUB_END, 0.06, Z_TH), (0, 0, 0), grey_d)
    sheave = rc.cylinder("moving_sheave", SHEAVE_R, Z_TH, (0, 0, 0), grey_m)

    # fixed tackle mount (anchor bracket) + fixed redirect sheave + anchor dot
    rc.box("tackle_bar", (0.10, 2.2, Z_TH), (ANCHOR[0], 5.55, 0), grey_d)
    rc.box("tackle_foot", (0.9, 0.10, Z_TH), (ANCHOR[0] - 0.45, 4.45, 0), grey_d)
    redirect = rc.cylinder("fixed_pulley", SHEAVE_R, Z_TH, (REDIRECT[0], REDIRECT[1], 0),
                           grey_m)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.07,
                                         location=(ANCHOR[0], ANCHOR[1], Z_TH / 2))
    anchor_dot = bpy.context.active_object
    anchor_dot.name = "anchor_dot"
    rc.set_material(anchor_dot, rc.make_material("anchor_m", rc.CABLE_BLUE))

    # insertion diamond rides the moving output tendon end
    ins = bpy.data.meshes.new("insertion")
    s = 0.11
    ins.from_pydata([(0, s, 0), (s, 0, 0), (0, -s, 0), (-s, 0, 0)], [],
                    [(0, 1, 2, 3)])
    ins.update()
    ins_ob = bpy.data.objects.new("insertion", ins)
    bpy.context.scene.collection.objects.link(ins_ob)
    rc.set_material(ins_ob, rc.make_material("ins_m", rc.CABLE_BLUE))
    ins_ob.location = (INSERTION0[0], INSERTION0[1], Z_TH / 2)

    # cable (keyframed polyline, 16 points)
    cable = rc.poly_curve("cable", cable_points(0.0), 0.032, material=cable_m)

    frames = list(range(0, N_FRAMES))
    deltas = [delta_of_frame(k) for k in frames]
    knees = [knee_pose(d) for d in deltas]
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.07, location=(*knees[0]['I'],Z_TH))
    pivot_marker=bpy.context.active_object
    rc.set_material(pivot_marker,rc.make_material('pivot_ink',rc.ARM_GOLD))
    rc.keyframe_linear(pivot_marker,'location',frames,[(*p['I'],Z_TH) for p in knees])
    leader=rc.poly_curve('pivot_leader',[(4.8,2.9,LABEL_Z),(*knees[0]['I'],Z_TH)],0.015,material=grey_d)
    rc.keyframe_polyline(leader,frames,[[(4.8,2.9,LABEL_Z),(*p['I'],Z_TH)] for p in knees])
    tibia = bpy.data.objects['tibia']
    rc.keyframe_linear(tibia, 'location', frames, [(*p['t'], 0) for p in knees])
    rc.keyframe_linear(tibia, 'rotation_euler', frames, [p['theta'] for p in knees], index=2)
    for j, (fixed, moving) in enumerate(((pk.A,'D'),(pk.B,'C'))):
        rc.keyframe_polyline(bpy.data.objects[f'kneelink{j}'],frames,
                            [[(*fixed,0),(*p[moving],0)] for p in knees])
    for j, key in ((2,'C'),(3,'D')):
        rc.keyframe_linear(bpy.data.objects[f'kneedot{j}'],'location',frames,
                           [(*p[key],Z_TH/2) for p in knees])
    polys = [cable_points(d) for d in deltas]
    rc.keyframe_polyline(cable, frames, polys)

    # BPA bodies shrink from the right; crosshead + sheave shift left by delta
    for ob, yy in zip(bpas, BPA_Y):
        locs = [((BPA_X0 + (BPA_X1 - d)) / 2.0, yy, 0) for d in deltas]
        rc.keyframe_linear(ob, "location", frames, locs)
        for i in range(3):
            rc.keyframe_linear(ob, "scale", frames,
                               [((bpa_len0 - d) / bpa_len0, 1.0, 1.0)[i] for d in deltas],
                               index=i)
        rc.keyframe_linear(bpy.data.objects[f"stub{yy}"], "location", frames,
                           [(BPA_X1 - d + 0.14, yy, 0) for d in deltas])
    rc.keyframe_linear(cross, "location", frames,
                       [(STUB_END - d, sum(BPA_Y) / 2.0, 0) for d in deltas])
    rc.keyframe_linear(cross_arm, "location", frames,
                       [((SHEAVE0[0] - d + STUB_END - d) / 2.0, SHEAVE0[1], 0) for d in deltas])
    rc.keyframe_linear(sheave, "location", frames,
                       [(SHEAVE0[0] - d, SHEAVE0[1], 0) for d in deltas])
    rc.keyframe_linear(ins_ob, "location", frames,
                       [(*p['P'],Z_TH/2) for p in knees])

    # ---- motion arrows -----------------------------------------------------
    red_a = rc.make_material("arrow_red", rc.BPA_RED)
    blue_a = rc.make_material("arrow_blue", rc.CABLE_BLUE)
    bpa_arrows = []
    for yy in BPA_Y:
        shaft, head = rc.arrow(f"bpa_arrow{yy}", (2.3, yy - 0.55, LABEL_Z),
                               (1.6, yy - 0.55, LABEL_Z), 0.045, red_a)
        bpa_arrows.append((shaft, head, yy))
    out_shaft, out_head = rc.arrow("out_arrow", (TENDON_X + 0.55, 2.0, LABEL_Z),
                                   (TENDON_X + 0.55, 1.2, LABEL_Z), 0.045, blue_a)

    def update_arrows(k):
        d = deltas[k]
        for shaft, head, yy in bpa_arrows:
            tail = (2.3 - d, yy - 0.55, LABEL_Z)
            headv = (1.6 - d, yy - 0.55, LABEL_Z)
            rc.aim_arrow(shaft, head, tail, headv)
        p=knees[k]
        direction=(p['tangent']-p['P'])/np.linalg.norm(p['tangent']-p['P'])
        start=p['P']+np.array([0.50,-0.25])
        tail=(*start,LABEL_Z)
        headv=(*(start+0.8*direction),LABEL_Z)
        rc.aim_arrow(out_shaft, out_head, tail, headv)

    # ---- labels (flat, facing +z) -------------------------------------------
    def lab(body, name, size=0.24, color=rc.INK, align="CENTER"):
        ob = rc.text_label(body, name, size, color=color, align=align)
        ob.location = (ob.location.x, ob.location.y, LABEL_Z)
        return ob

    L = {}
    L["mount"] = lab("Origin mount", "L_mount", 0.22)
    L["mount"].location = (MOUNT_X, 7.35, LABEL_Z)
    L["bpas"] = lab("Two parallel BPAs", "L_bpas", 0.24, color=rc.BPA_RED)
    L["bpas"].location = (2.4, 7.15, LABEL_Z)
    L["moving"] = lab("Moving pulley", "L_moving", 0.22)
    L["moving"].location = (SHEAVE0[0] - 0.2, 6.95, LABEL_Z)
    L["fixedpul"] = lab("Fixed pulley", "L_fixedpul", 0.22)
    L["fixedpul"].location = (REDIRECT[0] + 1.15, 4.62, LABEL_Z)
    L["tackle"] = lab("Fixed tackle mount (femur)", "L_tackle", 0.20)
    L["tackle"].location = (ANCHOR[0] + 0.35, 6.9, LABEL_Z)
    L["tendon"] = lab("Output tendon", "L_tendon", 0.24, color=rc.CABLE_BLUE,
                      align="LEFT")
    L["tendon"].location = (TENDON_X + 0.45, 3.6, LABEL_Z)
    L["ins"] = lab("Tendon attachment", "L_ins", 0.20, align="RIGHT")
    rc.keyframe_linear(L["ins"], "location", frames,
                       [(p['P'][0] - 0.20, p['P'][1] - 0.30, LABEL_Z) for p in knees])
    L["femur"] = lab("Femur", "L_femur", 0.22, color=(0.35, 0.34, 0.32))
    L["femur"].location = (3.2, 4.85, LABEL_Z)
    L["tibia"] = lab("Tibia", "L_tibia", 0.22, color=(0.35, 0.34, 0.32))
    rc.keyframe_linear(L["tibia"], 'location', frames,
                       [(*(p['Q']@np.array([8.35,1.15])+p['t']),LABEL_Z) for p in knees])
    L["knee"] = lab("Moving pivot: link intersection", "L_knee", 0.19, color=(0.35, 0.34, 0.32),
                    align="RIGHT")
    L["knee"].location = (5.9, 2.75, LABEL_Z)
    L["force"] = lab("Ideal: output force = half the\ncombined BPA force",
                     "L_force", 0.21, color=rc.BPA_RED, align="LEFT")
    L["force"].location = (0.55, 1.5, LABEL_Z)
    L["contr"] = lab("BPA contraction (each)", "L_contr", 0.19, color=rc.BPA_RED,
                     align="CENTER")
    L["contr"].location = (1.95, 4.78, LABEL_Z)
    L["travel"] = lab("Output travel = 2 x\nBPA contraction", "L_travel", 0.17,
                      color=rc.CABLE_BLUE, align="LEFT")
    L["travel"].location = (0.55, 0.8, LABEL_Z)

    cam_data = bpy.data.cameras.new("cam")
    cam_data.type = "ORTHO"
    cam_data.ortho_scale = 12.0
    cam = bpy.data.objects.new("camera", cam_data)
    bpy.context.scene.collection.objects.link(cam)
    cam.location = (5.0, 4.0, 12.0)
    cam.rotation_euler = rc.look_rotation(cam.location, (5.0, 4.0, 0.0))
    scene.camera = cam

    rc.save_blend(rc.OUTPUT / "scenes" / "visual2_reverse_pulley.blend")

    if args.stills_only:
        k = int(np.argmax(deltas))
        scene.render.resolution_x = 1200
        scene.render.resolution_y = 800
        update_arrows(k)
        scene.frame_set(k)
        rc.render_frame(scene, k, str(rc.OUTPUT / "visual2" / "still_labeled.png"))
        # unlabeled copy: hide all text objects
        texts = [o for o in bpy.data.objects if o.type == "FONT"]
        for o in texts:
            o.hide_render = True
        for a in bpa_arrows:
            a[0].hide_render = True
            a[1].hide_render = True
        out_shaft.hide_render = True
        out_head.hide_render = True
        rc.render_frame(scene, k, str(rc.OUTPUT / "visual2" / "still_unlabeled.png"))
        print("STILLS_DONE")
        return

    outdir = rc.OUTPUT / "visual2" / "frames"
    outdir.mkdir(parents=True, exist_ok=True)
    for k in frames:
        update_arrows(k)
        rc.render_frame(scene, k, str(outdir / f"f{k:04d}.png"))
        if k % 30 == 0:
            print(f"rendered {k}/{N_FRAMES}")
    print("VISUAL2_DONE")


if __name__ == "__main__":
    main()
