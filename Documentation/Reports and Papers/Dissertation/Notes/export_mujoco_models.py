"""Render existing MuJoCo models for preliminary Methods figures.

Run with conda run --live-stream -n myoconv python <this file>.
No simulation, conversion, control callback, or model-file changes occur.
The converted robot uses saved keyframe 0. The synthetic knee is posed at
-35 degrees solely to expose the two actual tendon routes.
"""
from pathlib import Path
import hashlib
import json
import sys
import time
import xml.etree.ElementTree as ET

import mujoco
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parents[1] / "ProofFinal/figs/Preliminary"
sys.path.insert(0, str(ROOT / "Code/MuJoCo_SNS"))
from sns_bpa_demo import MJCF


def render(model, data, name, lookat, distance, azimuth, elevation):
    model.vis.global_.offwidth = 1600
    model.vis.global_.offheight = 1500
    model.vis.headlight.ambient[:] = [0.18, 0.18, 0.18]
    model.vis.headlight.diffuse[:] = [0.5, 0.5, 0.5]
    model.vis.headlight.specular[:] = [0.15, 0.15, 0.15]
    camera = mujoco.MjvCamera()
    camera.lookat[:] = lookat
    camera.distance = distance
    camera.azimuth = azimuth
    camera.elevation = elevation
    options = mujoco.MjvOption()
    options.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = True
    options.flags[mujoco.mjtVisFlag.mjVIS_CONSTRAINT] = False
    options.sitegroup[:] = 0
    # MuJoCo 2.3.7's Renderer predates close()/context-manager support.
    renderer = mujoco.Renderer(model, height=1500, width=1000)
    renderer.update_scene(data, camera=camera, scene_option=options)
    renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = True
    Image.fromarray(renderer.render()).save(OUT / name)
    return dict(lookat=lookat, distance=distance, azimuth=azimuth, elevation=elevation)


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    robot_path = ROOT / "Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_robot/gait2392_robot_cvt3.xml"
    robot = mujoco.MjModel.from_xml_path(str(robot_path))
    rd = mujoco.MjData(robot)
    mujoco.mj_resetDataKeyframe(robot, rd, 0)
    mujoco.mj_forward(robot, rd)
    # Display-only appearance settings, preserving geometry and tendon paths.
    robot.tendon_width[:] = 0.003
    robot.tendon_rgba[:] = [0.78, 0.12, 0.11, 1]
    robot.geom_rgba[robot.geom("ground-plane").id] = [0.92, 0.93, 0.95, 1]
    if "--viewer" in sys.argv:
        # Explicit opt-in: open the actual native GUI for documentary capture.
        # This loop synchronizes display only; it never calls mj_step.
        from mujoco import viewer as viewer_api
        with viewer_api.launch_passive(robot, rd) as viewer:
            viewer.cam.lookat[:] = [0, 0, 0.9]
            viewer.cam.distance = 2.8
            viewer.cam.azimuth = 215
            viewer.cam.elevation = -8
            started = time.monotonic()
            while viewer.is_running() and time.monotonic() - started < 300:
                viewer.sync()
                time.sleep(0.05)
        print("Viewer closed; no simulation steps executed.")
        return
    robot_camera = render(robot, rd, "mujoco_converted_robot.png", [0, 0, 0.90], 2.8, 215, -8)

    display_xml = ET.fromstring(MJCF)
    asset = ET.SubElement(display_xml, "asset")
    ET.SubElement(asset, "texture", type="skybox", builtin="gradient", rgb1="0.96 0.97 0.98", rgb2="0.96 0.97 0.98", width="256", height="1536")
    knee = mujoco.MjModel.from_xml_string(ET.tostring(display_xml, encoding="unicode"))
    kd = mujoco.MjData(knee)
    kd.qpos[knee.jnt_qposadr[knee.joint("knee").id]] = np.deg2rad(-35)
    knee.tendon_width[:] = 0.004
    knee.tendon_rgba[0] = [0.05, 0.55, 0.70, 1]
    knee.tendon_rgba[1] = [0.96, 0.65, 0.12, 1]
    knee.geom_rgba[0] = [0.48, 0.53, 0.6, 0.35]
    knee.geom_rgba[1] = [0.7, 0.34, 0.34, 0.65]
    mujoco.mj_forward(knee, kd)
    knee_camera = render(knee, kd, "mujoco_synthetic_knee.png", [-0.025, 0, 0.27], 0.9, 90, -5)
    info = dict(
        mujoco_version=mujoco.__version__,
        robot=dict(source=str(robot_path.relative_to(ROOT)), sha256=hashlib.sha256(robot_path.read_bytes()).hexdigest(),
                   keyframe=0, nbody=robot.nbody, nq=robot.nq, ntendon=robot.ntendon, nu=robot.nu, camera=robot_camera,
                   interpretation="Existing converted robot model, no SNS control or dynamics demonstrated."),
        knee=dict(source="Code/MuJoCo_SNS/sns_bpa_demo.py:MJCF", mjcf_sha256=hashlib.sha256(MJCF.encode()).hexdigest(),
                  angle_degrees=-35, nq=knee.nq, ntendon=knee.ntendon, nu=knee.nu, camera=knee_camera,
                  tendon_lengths_m=kd.ten_length.tolist(), tendon_names=[knee.tendon(i).name for i in range(knee.ntendon)],
                  interpretation="Illustrative static pose, not simulation output. Cyan: named flexor tendon; gold: named extensor tendon."),
        modifications="In-memory render settings only: camera, lighting, ground color, tendon colors and visual widths; knee background skybox and partially transparent thigh/shank to expose routes. Source geometry/physics unchanged; no steps executed.")
    (OUT / "mujoco_render_provenance.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
