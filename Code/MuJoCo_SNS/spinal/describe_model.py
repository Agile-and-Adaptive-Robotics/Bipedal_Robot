"""Probe the converted gait2392_simbody MJCF: joints, actuators, keyframes.

Run from Code/MuJoCo_SNS with the myo env:
    D:/Anaconda/envs/myo/python.exe spinal/describe_model.py
"""
import sys
from pathlib import Path

import mujoco

MODEL = Path(r"D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody"
             r"\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml")


def main():
    model = mujoco.MjModel.from_xml_path(str(MODEL))
    print(f"model: {MODEL.name}")
    print(f"nq={model.nq} nv={model.nv} nu={model.nu} nbody={model.nbody} "
          f"nkey={model.nkey} timestep={model.opt.timestep}")

    print("\n-- joints (name, type, nq) --")
    for j in range(model.njnt):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, j)
        print(f"  {name:28s} type={model.jnt_type[j]} nq={model.jnt_qposadr[j+1]-model.jnt_qposadr[j] if j+1 < model.njnt else model.nq-model.jnt_qposadr[j]}")

    print("\n-- actuators (idx, name, ctrlrange) --")
    for a in range(model.nu):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, a)
        lo, hi = model.actuator_ctrlrange[a]
        print(f"  {a:3d} {name:24s} [{lo:.2f}, {hi:.2f}]")

    print("\n-- keyframes --")
    for k in range(model.nkey):
        print(f"  key {k}: name={mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_KEY, k)!r} "
              f"time={model.key_time[k]:.3f}")


if __name__ == "__main__":
    sys.exit(main())
