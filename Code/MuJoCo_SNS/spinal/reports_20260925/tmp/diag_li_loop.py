import io, os, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
os.environ.setdefault("CONDA_PREFIX", r"C:\Users\Ben Bolen\.conda\envs\myo")
import numpy as np
import mujoco
HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
sys.path.insert(0, HERE)
from build_li_net import build, V_REST

m = mujoco.MjModel.from_xml_path(os.path.join(HERE, "w2l_mjcf.xml"))
d = mujoco.MjData(m)
net = build(dt=0.0002)
act = {a: mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
       for a in net.muscle_outputs}
jnt = {j: m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, j)]
       for j in ("hip_L", "knee_L", "hip_R", "knee_R")}
gids = [mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_GEOM, g) for g in
        ("foot_L_contact", "toe_L_contact", "foot_R_contact", "toe_R_contact")]

def watch_cells():
    names = ("L_CPG stance", "R_CPG stance", "L_hip stance PF", "L_hip swing PF",
             "R_hip stance PF", "R_hip swing PF", "L_hip stance MN MV",
             "L_hip swing MN MV")
    return [net.idx[n] for n in names], names

idxs, names = watch_cells()
u = net.make_inputs()
for i in range(3000):                    # 3 s
    t = i * 0.001
    c = np.zeros(4)
    for k in range(d.ncon):
        con = d.contact[k]
        for kk, gid in enumerate(gids):
            if con.geom1 == gid or con.geom2 == gid:
                f6 = np.zeros(6); mujoco.mj_contactForce(m, d, k, f6)
                if f6[0] > 1e-6: c[kk] = 1.0
    u[:] = 0.0
    net.set_all_tonics(u, extra={"R tonic drive": 3.0} if t < 0.2 else None)
    net.set_port(u, "L_foot ground contact", net.contact_current(None, c[0]))
    net.set_port(u, "L_toe ground contact", net.contact_current(None, c[1]))
    net.set_port(u, "R_foot ground contact", net.contact_current(None, c[2]))
    net.set_port(u, "R_toe ground contact", net.contact_current(None, c[3]))
    net.set_port(u, "L_hip middle", net.hip_current(None, np.degrees(d.qpos[jnt["hip_L"]])))
    net.set_port(u, "R_hip middle", net.hip_current(None, np.degrees(d.qpos[jnt["hip_R"]])))
    for _ in range(5):
        V = net.step(u)
    ctrl = net.muscle_ctrl(V)
    for a, v in ctrl.items():
        d.ctrl[act[a]] = v
    mujoco.mj_step(m, d)
    if i % 250 == 0:
        vs = " ".join(f"{V[k]:6.1f}" for k in idxs)
        ct = " ".join(f"{net.muscle_ctrl(V)[a]:.2f}" for a in
                      ("hip_L_ext", "hip_L_flx", "knee_L_ext", "ankle_L_ext"))
        print(f"t={t:4.1f} z={d.qpos[2]:5.2f} c={c.astype(int)} hipL={np.degrees(d.qpos[jnt['hip_L']]):6.1f} "
              f"kneeL={np.degrees(d.qpos[jnt['knee_L']]):6.1f} | V: {vs} | ctrl: {ct}")
print("cells:", names)
