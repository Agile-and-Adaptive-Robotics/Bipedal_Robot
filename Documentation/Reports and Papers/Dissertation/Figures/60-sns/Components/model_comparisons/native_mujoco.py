import argparse
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("CONDA_PREFIX", str(Path(sys.executable).parent))
import mujoco
import numpy as np
from PIL import Image


def render_native(model_path, output, azimuth, elevation, distance, height, bridge=None):
    model = mujoco.MjModel.from_xml_path(str(model_path))
    state = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, state, 0)
    if bridge:
        expected = "4f998681d610e20073fc077f9c477d7f720e54061e9a594c241f34d9b65c8870"
        if hashlib.sha256(model_path.read_bytes()).hexdigest() != expected:
            raise ValueError("The recorded robotbody and supplied model differ.")
        with gzip.open(bridge, "rt") as stream:
            record = json.load(stream)
        selected = model.actuator(record["actuator"]).id
        addresses = [model.joint(name).qposadr[0] for name in record["jointNames"]]
        for condition in ("off", "on"):
            samples = record[condition]
            trial = mujoco.MjData(model)
            mujoco.mj_resetDataKeyframe(model, trial, 0)
            for k in range(1, len(samples["t"])):
                trial.ctrl[selected] = samples["command"][k-1]
                mujoco.mj_step(model, trial)
                np.testing.assert_allclose(trial.actuator_force[selected], samples["force"][k], atol=1e-10, rtol=0)
            for k in range(1, len(samples["t"])):
                trial.qpos[addresses] = samples["q"][k]
                mujoco.mj_forward(model, trial)
                np.testing.assert_allclose(trial.actuator_length[selected], samples["length"][k], atol=1e-10, rtol=0)
        frame, = np.flatnonzero(np.isclose(record["on"]["t"], .025))
        state.qpos[addresses] = record["on"]["q"][frame]
    mujoco.mj_forward(model, state)
    camera = mujoco.MjvCamera()
    mujoco.mjv_defaultCamera(camera)
    camera.type = mujoco.mjtCamera.mjCAMERA_FREE
    camera.azimuth = azimuth
    camera.elevation = elevation
    camera.distance = distance
    camera.lookat[:] = [-0.055, 0, height]
    width, pixels = 720, 1080
    model.vis.global_.offwidth = width
    model.vis.global_.offheight = pixels
    renderer = mujoco.Renderer(model, pixels, width)
    renderer.update_scene(state, camera)
    Image.fromarray(renderer.render()).save(output)
    print('Angles and translation',[(j, state.qpos[model.joint(j).qposadr[0]]) for j in ['pelvis_ty','hip_flexion_r','knee_angle_r','ankle_angle_r']])
    print('tendon colors',np.unique(model.tendon_rgba,axis=0),'width',np.unique(model.tendon_width))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('model',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--azimuth',type=float,default=155)
    p.add_argument('--elevation',type=float,default=-6)
    p.add_argument('--distance',type=float,default=2.36)
    p.add_argument('--height',type=float,default=.85)
    p.add_argument('--bridge',type=Path)
    args = p.parse_args()
    render_native(args.model,args.output,args.azimuth,args.elevation,args.distance,args.height,args.bridge)
