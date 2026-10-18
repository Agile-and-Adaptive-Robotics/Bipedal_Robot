"""Isometric print render of the converted Gait2392 robot-body model.

Ask (2026-10-18): isometric view (elevated ~30 deg, azimuth ~45 deg) of the
converted robot-body musculoskeletal model at its saved reference pose,
muscle (tendon) paths visible, white background, 300 dpi PNG, >=1800 px
long edge. Display-only in-memory settings; no physics stepped, no file
changes. Pattern follows Documentation/Reports and Papers/Dissertation/
Notes/export_mujoco_models.py (keyframe 0, mjVIS_TENDON, Renderer).
"""
import json
import sys
from pathlib import Path

import mujoco
import numpy as np
from PIL import Image

ROOT = Path(r"D:\Github\Bipedal_Robot")
XML = ROOT / r"Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_robot\gait2392_robot_cvt3.xml"
OUT_PNG = ROOT / r"Documentation\Reports and Papers\Dissertation\Figures\30-torque\opensim_robot_isometric.png"

W, H = 1440, 2160          # final long edge 2160 px >= 1800; 2:3 portrait
SS = 2                     # supersample 2x (2880x4320), LANCZOS downsample for print-crisp AA
DPI = (300, 300)
AZIMUTH = 45.0             # ask: azimuth ~45 deg
ELEVATION = -30.0          # MuJoCo: negative elevation = camera ABOVE horizon -> 30 deg elevated view
LOOKAT = [0.0, 0.0, 0.90]
DISTANCE = 2.8


def main():
    model = mujoco.MjModel.from_xml_path(str(XML))
    data = mujoco.MjData(model)
    mujoco.mj_resetDataKeyframe(model, data, 0)   # saved reference pose "default-pose"
    mujoco.mj_forward(model, data)

    # ---- display-only, in-memory appearance settings (no model-file edits) ----
    model.tendon_width[:] = 0.003
    model.tendon_rgba[:] = [0.78, 0.12, 0.11, 1]              # muscle paths in deep red
    # White background (A/B-probed in ab_background.py): the model ships a gradient
    # SKYBOX texture (xml line 23) that tints the sky; whitening tex_rgb renders it
    # BLACK and disabling mjRND_SKYBOX exposes the black clear color. Instead, enlarge
    # the ground plane in memory to 40x40 m so the 30-deg-down camera never sees past
    # it, and color the plane pure white -> every background pixel is white.
    ground_id = model.geom("ground-plane").id
    model.geom_size[ground_id] = [40, 40, 0.125]
    model.geom_rgba[ground_id] = [1, 1, 1, 1]
    model.vis.rgba.fog[:] = [1, 1, 1, 1]
    model.vis.global_.offwidth = W * SS                       # enlarge offscreen framebuffer
    model.vis.global_.offheight = H * SS
    model.vis.headlight.ambient[:] = [0.18, 0.18, 0.18]
    model.vis.headlight.diffuse[:] = [0.5, 0.5, 0.5]
    model.vis.headlight.specular[:] = [0.15, 0.15, 0.15]

    camera = mujoco.MjvCamera()
    camera.lookat[:] = LOOKAT
    camera.distance = DISTANCE
    camera.azimuth = AZIMUTH
    camera.elevation = ELEVATION

    options = mujoco.MjvOption()
    options.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = True       # muscle paths visible
    options.flags[mujoco.mjtVisFlag.mjVIS_CONSTRAINT] = False
    options.sitegroup[:] = 0                                   # hide path-point sites

    renderer = mujoco.Renderer(model, height=H * SS, width=W * SS)
    renderer.update_scene(data, camera=camera, scene_option=options)
    renderer.scene.flags[mujoco.mjtRndFlag.mjRND_SHADOW] = True
    arr = renderer.render()

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    img = Image.fromarray(arr)
    if SS > 1:
        img = img.resize((W, H), Image.LANCZOS)
    img.save(OUT_PNG, dpi=DPI)
    print("saved:", OUT_PNG)
    print(json.dumps(dict(
        mujoco_version=mujoco.__version__,
        model=str(XML.relative_to(ROOT)),
        keyframe=mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_KEY, 0),
        nbody=int(model.nbody), nq=int(model.nq),
        ntendon=int(model.ntendon), nu=int(model.nu),
        camera=dict(lookat=LOOKAT, distance=DISTANCE,
                    azimuth=AZIMUTH, elevation=ELEVATION),
        pixels=[W, H], dpi=list(DPI), supersample=SS,
        ground_plane_enlarged_display_only=True,
        engine="mujoco.Renderer offscreen (native WGL context)",
    ), indent=2))


if __name__ == "__main__":
    sys.exit(main())
