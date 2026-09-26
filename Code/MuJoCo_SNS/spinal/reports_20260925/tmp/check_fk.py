# -*- coding: utf-8 -*-
"""FK sanity: key anchor world positions, muscle rest distances, tibia convention."""
import sys, os
sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco")
import make_w2l_mjcf as M

m = M.parse_aproj()
B = {b.name: b for b in m["bodies"]}

print("== key world positions (AL frame, y-up, meters) ==")
for nm in ("Root", "femur_L", "tibia_L", "foot_L", "toe_L", "femur_R", "tibia_R", "foot_R", "toe_R"):
    b = [x for x in m["bodies"] if x.name == nm and x.type == "Box"][0]
    print("  %-8s pos=(%.4f, %.4f, %.4f)" % (nm, *b.p_world))

print("\n== joint anchors ==")
for b in m["bodies"]:
    if b.joint:
        j = b.joint
        anchor = M.vadd(b.p_world, M.mat_vec(b.R_world, j["pos_local_al"]))
        axis = M.mat_vec(M.mat_mul(b.R_world, M.mat_mul(b.R_al, j["R_local_al"])), [1, 0, 0])
        print("  %-8s anchor=(%.4f, %.4f, %.4f) axis=(%.3f, %.3f, %.3f) limits=[%s, %s] deg" %
              (j["name"], *anchor, *axis, j["lower_deg"], j["upper_deg"]))

print("\n== muscle attachment distances vs RestingLength ==")
for mus in m["muscles"]:
    pts = []
    for aid in mus["attach_ids"]:
        att = m["by_id"][aid]
        pts.append(att.site_pos_world_al)
    L = sum(M.vnorm(M.vsub(pts[i + 1], pts[i])) for i in range(len(pts) - 1))
    straight = M.vnorm(M.vsub(pts[-1], pts[0]))
    print("  %-12s L0=%.3f Lw=%.3f path=%.3f straight=%.3f attach=%s" %
          (mus["name"], mus["resting_length"], mus["lwidth"], L, straight,
           [m["by_id"][a].name for a in mus["attach_ids"]]))

print("\n== foot contact boxes ==")
for nm in ("foot_L_contact", "toe_L_contact", "foot_R_contact", "toe_R_contact"):
    b = [x for x in m["bodies"] if x.name == nm][0]
    print("  %-16s world=(%.4f, %.4f, %.4f) lwh=%s parent=%s" %
          (nm, *b.p_world, b.lwh, b.parent.name))

print("\n== tibia_L convention check ==")
t = [x for x in m["bodies"] if x.name == "tibia_L" and x.type == "Box"][0]
p_lr, r_lr = M.parse_pos_rot_rb = None, None
import xml.etree.ElementTree as ET
tree = ET.parse(M.APROJ)
org = tree.getroot().find(".//Organism")
for rb in org.iter("RigidBody"):
    if M.txt(rb, "Name") == "tibia_L" and M.txt(rb, "Type") == "Box":
        R_lm, p_lm = M.parse_local_matrix(rb)
        p_lr, r_lr = M.parse_pos_rot(rb)
        R_er = M.euler_xyz_deg(*r_lr)
        print("  LocalMatrix R:", [[round(c, 5) for c in row] for row in R_lm])
        print("  Euler RzRyRx :", [[round(c, 5) for c in row] for row in R_er])
        print("  rot=(%s) dev=%.4g" % (r_lr, max(abs(R_lm[i][j] - R_er[i][j]) for i in range(3) for j in range(3))))
        # try other orders
        for nm, fn in (("RyRx", lambda: M.mat_mul(M.ry(r_lr[1] * M.D2R), M.rx(r_lr[0] * M.D2R))),
                       ("RxRyRz", lambda: M.mat_mul(M.rx(r_lr[0] * M.D2R), M.mat_mul(M.ry(r_lr[1] * M.D2R), M.rz(r_lr[2] * M.D2R)))),
                       ("RzRxRy", lambda: M.mat_mul(M.rz(r_lr[2] * M.D2R), M.mat_mul(M.rx(r_lr[0] * M.D2R), M.ry(r_lr[1] * M.D2R))))):
            R = fn()
            dev = max(abs(R_lm[i][j] - R[i][j]) for i in range(3) for j in range(3))
            print("    order %-7s dev=%.4g" % (nm, dev))
        break

print("\n== spring endpoints ==")
for spr in m["springs"]:
    for aid in spr["attach_ids"]:
        att = m["by_id"][aid]
        print("  %s: %s at (%.4f, %.4f, %.4f) parent=%s" %
              (spr["name"], att.name, *att.site_pos_world_al, att.parent.name))
