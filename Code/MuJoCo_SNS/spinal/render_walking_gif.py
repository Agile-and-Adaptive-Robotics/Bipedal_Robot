"""Render skeletal walking gifs from spinal run npz using MuJoCo offscreen renderer.
v3: handles the gait2392_simbody model structure (3 slide + 3 hinge pelvis joints,
no free joint), drives pelvis from COM data, and maps all 13 recorded joints."""
import sys
from pathlib import Path

import numpy as np

SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")


def render_walking_gif(npz_path, output_gif, model_xml=None, fps=30, seconds=6):
    import mujoco
    from PIL import Image

    if model_xml is None:
        model_xml = (SPINAL.parent.parent / "Solid_Models" / "OpenSim" /
                     "Gait2392_Robotbody" / "mjc" / "gait2392_simbody" /
                     "gait2392_simbody_cvt3.xml")
    model_xml = Path(model_xml)
    if not model_xml.exists():
        print(f"ERROR: model xml not found: {model_xml}")
        return False

    d = np.load(npz_path, allow_pickle=True)
    q = d["q"]
    t = d["t"] if "t" in d else np.arange(len(q)) * 0.002
    com = d["com"] if "com" in d else None
    key_joints = [str(x) for x in d["key_joints"]] if "key_joints" in d else []

    model = mujoco.MjModel.from_xml_path(str(model_xml))
    data = mujoco.MjData(model)

    # Build the joint name -> qpos address map
    name2adr = {}
    for i in range(model.njnt):
        name = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, i)
        name2adr[name] = model.jnt_qposadr[i]

    # The pelvis slide joints (tx, ty, tz) are the FIRST 3 qpos
    # The pelvis hinge joints (tilt, list, rotation) are qpos 3-5
    # These are in key_joints as well (the runner records them)
    # Map ALL recorded joints
    joint_map = {}
    for i, rec_name in enumerate(key_joints):
        if rec_name in name2adr:
            joint_map[i] = (rec_name, name2adr[rec_name])

    print(f"Model: {model.nq} qpos, {model.ngeom} geoms")
    print(f"Mapped {len(joint_map)}/{len(key_joints)} joints")
    print(f"Pelvis tx={name2adr.get('pelvis_tx')}, ty={name2adr.get('pelvis_ty')}, tz={name2adr.get('pelvis_tz')}")

    # Setup renderer
    renderer = mujoco.Renderer(model, height=480, width=640)
    camera = mujoco.MjvCamera()
    camera.lookat[:] = [0, 0, 0.85]
    camera.distance = 3.5
    camera.elevation = -8
    camera.azimuth = 90

    # Frames
    total = len(q)
    target = int(seconds * fps)
    step = max(1, total // target)
    indices = list(range(0, total, step))[:target]
    skip = len(indices) // 10
    indices = indices[skip:]

    # Pelvis slide qpos addresses
    tx_adr = name2adr.get("pelvis_tx", 0)
    ty_adr = name2adr.get("pelvis_ty", 1)
    tz_adr = name2adr.get("pelvis_tz", 2)

    images = []
    for fi in indices:
        mujoco.mj_resetData(model, data)

        # Set pelvis translations from COM.
        # The converter preserved OpenSim joint names but the MuJoCo world frame
        # is z-up. Empirically: pelvis_ty (qpos[1]) moves the body in world Z
        # (up), pelvis_tz (qpos[2]) moves it in world -Y (lateral). So:
        #   forward  = com[:,0] -> qpos[0] (pelvis_tx)
        #   up       = com[:,2] -> qpos[1] (pelvis_ty, which is MuJoCo Z-up)
        #   lateral  = com[:,1] -> qpos[2] (pelvis_tz, which is MuJoCo -Y)
        if com is not None and fi < len(com):
            data.qpos[tx_adr] = com[fi, 0]  # forward
            data.qpos[ty_adr] = com[fi, 2]  # UP (world Z, from COM column 2)
            data.qpos[tz_adr] = -com[fi, 1]  # lateral (world -Y, negated COM column 1)

        # Set all recorded hinge joints
        for col, (name, adr) in joint_map.items():
            if col < q.shape[1]:
                data.qpos[adr] = np.radians(q[fi, col])

        mujoco.mj_forward(model, data)

        # Track the model with the camera (follow pelvis height)
        camera.lookat[2] = data.qpos[tz_adr] * 0.5 + 0.4

        renderer.update_scene(data, camera=camera)
        img = renderer.render()
        images.append(Image.fromarray(img))

    if not images:
        return False

    out = Path(output_gif)
    out.parent.mkdir(parents=True, exist_ok=True)
    images[0].save(out, save_all=True, append_images=images[1:], duration=1000 // fps, loop=0)
    print(f"Saved {len(images)} frames to {out}")
    return True


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        ok = render_walking_gif(sys.argv[1], sys.argv[2],
                                sys.argv[3] if len(sys.argv) > 3 else None)
        sys.exit(0 if ok else 1)
    print("Usage: python render_walking_gif.py <npz> <output.gif> [model.xml]")
