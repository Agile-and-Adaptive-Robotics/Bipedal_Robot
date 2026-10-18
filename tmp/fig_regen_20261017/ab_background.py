"""A/B experiment: find what paints the dark top band, small renders."""
import numpy as np
import mujoco
from PIL import Image

XML = r"D:\Github\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_robot\gait2392_robot_cvt3.xml"

flags = [f for f in dir(mujoco.mjtRndFlag) if f.startswith("mjRND")]
print("mjtRndFlag members:", flags)
for f in flags:
    print(" ", f, "=", int(getattr(mujoco.mjtRndFlag, f)))


def probe(tag, skybox=None, fog=None, haze=None, whiten=False, plane_big=False):
    m = mujoco.MjModel.from_xml_path(XML)
    d = mujoco.MjData(m)
    mujoco.mj_resetDataKeyframe(m, d, 0)
    mujoco.mj_forward(m, d)
    m.vis.global_.offwidth = 320
    m.vis.global_.offheight = 240
    m.vis.rgba.fog[:] = [1, 1, 1, 1]
    m.geom_rgba[m.geom("ground-plane").id] = [0.94, 0.95, 0.97, 1]
    if whiten:
        n = int(m.tex_height[0]) * int(m.tex_width[0]) * 3
        m.tex_rgb[0:n] = 1.0
    if plane_big:
        m.geom_size[m.geom("ground-plane").id] = [40, 40, 0.125]
    cam = mujoco.MjvCamera()
    cam.lookat[:] = [0, 0, 0.9]
    cam.distance = 2.8
    cam.azimuth = 45
    cam.elevation = -30
    opt = mujoco.MjvOption()
    opt.flags[mujoco.mjtVisFlag.mjVIS_TENDON] = True
    opt.sitegroup[:] = 0
    r = mujoco.Renderer(m, height=240, width=320)
    r.update_scene(d, camera=cam, scene_option=opt)
    if skybox is not None:
        r.scene.flags[mujoco.mjtRndFlag.mjRND_SKYBOX] = skybox
    if fog is not None:
        r.scene.flags[mujoco.mjtRndFlag.mjRND_FOG] = fog
    if haze is not None:
        r.scene.flags[mujoco.mjtRndFlag.mjRND_HAZE] = haze
    a = np.asarray(Image.fromarray(r.render()).convert("RGB")).astype(int)
    pts = {"tl": a[3, 3].tolist(), "tmid": a[3, 160].tolist(), "r30": a[30, 160].tolist(),
           "r60": a[60, 160].tolist(), "c": a[120, 160].tolist(), "bl": a[236, 316].tolist()}
    print(f"{tag:34s} {pts}")


probe("A fresh (fog-white only)")
probe("B whiten tex_rgb", whiten=True)
probe("C skybox flag OFF", skybox=False)
probe("D skybox OFF + whiten", skybox=False, whiten=True)
probe("E haze OFF", haze=False)
probe("F fog OFF", fog=False)
probe("G plane 40m", plane_big=True)
probe("H skybox OFF + plane 40m", skybox=False, plane_big=True)
