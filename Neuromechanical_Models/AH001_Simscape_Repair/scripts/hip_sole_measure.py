#!/usr/bin/env python3
"""hip_sole_measure.py -- neutral-pose hip->sole vertical distance for
(1) the converted OpenSim Gait2392 human model and (2) the AH001 robot model,
both loaded in MuJoCo. Vertical = the model's gravity axis. Ben's tibia
directive: shorten the tibia so the ROBOT's neutral hip->sole matches OpenSim's.
"""
import os
import sys

import mujoco  # noqa: E402

REPO = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot"
MJC_DIR = os.path.join(REPO, r"Solid_Models\OpenSim\Gait2392_Robotbody\mjc")
ROBOT = os.path.join(REPO, r"Code\MuJoCo_SNS\robot_torso\ah001_neural.xml")
OUT = os.path.join(REPO, r"Neuromechanical_Models\AH001_Simscape_Repair\logs\hip_sole_measure_20261102.json")


def geom_low_y(m, d, gid):
    """lowest point of a geom along the vertical (gravity) axis."""
    gtype = m.geom_type[gid]
    pos = d.geom_xpos[gid]
    mat = d.geom_xmat[gid].reshape(3, 3)
    size = m.geom_size[gid]
    if gtype == mujoco.mjtGeom.mjGEOM_SPHERE:
        return pos - size[0]
    if gtype == mujoco.mjtGeom.mjGEOM_BOX:
        # half-extents in geom frame
        ext = size[:3] * 0 + size[:3]
        # world extent along vertical axis (world up = -grav): project each axis
        low = pos.copy()
        for i in range(3):
            low -= abs(mat[:, i]) * ext[i]
        return low
    if gtype == mujoco.mjtGeom.mjGEOM_CYLINDER:
        low = pos.copy()
        for i in range(2):  # x, z axes of cylinder frame
            low -= abs(mat[:, i]) * size[i]
        low -= abs(mat[:, 1]) * size[1] * 0  # height axis handled below
        low[1] -= size[1]  # approx: half height along local y? conservative
        return low
    if gtype == mujoco.mjtGeom.mjGEOM_CAPSULE:
        # cylinder + spheres; axis = local z
        axis = mat[:, 2]
        tip = size[1] + size[0]
        low = pos.copy()
        low -= abs(axis) * size[1]
        low -= size[0]  # cap radius on all axes
        return low
    if gtype == mujoco.mjtGeom.mjGEOM_MESH:
        # conservative bounding radius
        return pos - m.geom_rbound[gid]
    return pos - m.geom_rbound[gid]


def measure(path, hip_body_pat, foot_body_pats, label):
    m = mujoco.MjModel.from_xml_path(path)
    d = mujoco.MjData(m)
    mujoco.mj_resetData(m, d)
    mujoco.mj_forward(m, d)
    grav = m.opt.gravity
    up = -grav / (grav @ grav) ** 0.5  # unit vertical (up)
    # hip joint anchors: all joints of the hip bodies
    hips = []
    for j in range(m.njnt):
        bid = m.jnt_bodyid[j]
        bname = mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_BODY, bid)
        if hip_body_pat(bname):
            hips.append((mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, j),
                         d.xanchor[j].copy()))
    # sole: lowest point over foot-body geoms
    soles = []
    for g in range(m.ngeom):
        bid = m.geom_bodyid[g]
        bname = mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_BODY, bid)
        if foot_body_pats(bname):
            low = geom_low_y(m, d, g)
            soles.append((bname, g, low))
    res = dict(label=label, up=up.tolist(), hips=[], sole_min=None)
    for nm, a in hips:
        res["hips"].append(dict(joint=nm, anchor=a.tolist(), h=float(a @ up)))
    best = None
    for bn, g, low in soles:
        h = float(low @ up)
        res.setdefault("sole_geoms", []).append(dict(body=bn, geom=g, h=h))
        if best is None or h < best:
            best = h
    res["sole_min"] = best
    res["hip_minus_sole"] = [h["h"] - best for h in res["hips"]]
    return res, m, d


if __name__ == "__main__":
    print("mjc dir:", os.listdir(MJC_DIR) if os.path.isdir(MJC_DIR) else "MISSING")
    results = {}
    # OpenSim human: hip = femur_r body joints; foot = calcn_r + toes_r
    xmls = [f for f in os.listdir(MJC_DIR) if f.endswith(".xml")] if os.path.isdir(MJC_DIR) else []
    print("xmls:", xmls)
    main_xml = None
    for cand in ("gait2392_robot.xml", "gait2392.xml", "converted_model.xml"):
        if cand in xmls:
            main_xml = cand
            break
    if main_xml is None and xmls:
        main_xml = xmls[0]
    if main_xml:
        r, _, _ = measure(os.path.join(MJC_DIR, main_xml),
                          lambda b: b in ("femur_r", "femur_l"),
                          lambda b: b in ("calcn_r", "calcn_l", "toes_r", "toes_l",
                                          "calcn_r_l", "foot_r", "foot_l"),
                          "opensim_gait2392")
        results["opensim"] = r
        print("OpenSim hip heights:", [h["h"] for h in r["hips"]], "sole:", r["sole_min"])
    # Robot
    r2, _, _ = measure(ROBOT,
                       lambda b: b in ("thigh_r", "thigh_l"),
                       lambda b: b in ("foot_r", "foot_l"),
                       "ah001_robot")
    results["robot"] = r2
    print("Robot hip heights:", [h["h"] for h in r2["hips"]], "sole:", r2["sole_min"])

    import json
    with open(OUT, "w") as f:
        json.dump(results, f, indent=1)
    # summary
    if "opensim" in results and results["opensim"]["hip_minus_sole"]:
        o = max(results["opensim"]["hip_minus_sole"])
        rb = max(results["robot"]["hip_minus_sole"])
        print(f"\nOpenSim hip->sole (neutral): {o:.4f} m")
        print(f"Robot   hip->sole (neutral): {rb:.4f} m")
        print(f"Tibia shortening needed: {rb - o:.4f} m ({(rb - o) * 1000:.1f} mm)")
