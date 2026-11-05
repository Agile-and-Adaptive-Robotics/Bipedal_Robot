"""render_common.py — shared Blender helpers for the defense teaching visuals.

World frame convention for ALL three scenes: the anatomical/OpenSim frame
(x anterior, y superior, z lateral-right). Blender's default Z-up is ignored;
the camera provides the view orientation. White background throughout.
"""
import math
import struct
from pathlib import Path
import numpy as np
import bpy
from mathutils import Vector, Matrix

# scripts/ -> teaching_visuals -> v21 -> Defense -> Dissertation -> Reports and
# Papers -> Documentation -> <repo root>
REPO_ROOT = Path(__file__).resolve().parents[7]
GEOM_DIR = REPO_ROOT / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody" / "mjc" / "gait2392_simbody" / "Geometry"
HERE = Path(__file__).resolve().parents[1]          # teaching_visuals_20261104
OUTPUT = HERE / "build_fixed_20261105"

# ---------------------------------------------------------------- palette
INK = (0.10, 0.10, 0.12)          # label ink (near black)
BONE = (0.55, 0.54, 0.52)         # bone gray
BONE_EDGE = (0.45, 0.44, 0.42)
BPA_RED = (0.55, 0.08, 0.20)
CABLE_BLUE = (0.00, 0.25, 0.50)
CAM_PURPLE = (0.45, 0.35, 0.60)   # shaped cam drum
PULLEY_GRAY = (0.55, 0.56, 0.58)  # pulleys
ARM_GOLD = (0.65, 0.35, 0.02)
AXIS_X = (0.85, 0.33, 0.00)       # Okabe-Ito
AXIS_Y = (0.00, 0.60, 0.45)
AXIS_Z = (0.00, 0.45, 0.70)


# ---------------------------------------------------------------- engine
def setup_scene(render_engine="auto", samples=24, res=(720, 720)):
    # All animated keys in these scenes interpolate LINEARLY (the motion is
    # computed per-frame, so easing is baked into the key values).
    try:
        bpy.context.preferences.edit.keyframe_new_interpolation_type = "LINEAR"
    except (AttributeError, TypeError):
        pass
    scene = bpy.context.scene
    # clear factory objects
    for ob in list(bpy.data.objects):
        bpy.data.objects.remove(ob, do_unlink=True)

    world = bpy.data.worlds.new("white")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes["Background"]
    bg.inputs[0].default_value = (1.0, 1.0, 1.0, 1.0)
    bg.inputs[1].default_value = 1.0    # strength 1.0 = true white background

    engine_used = render_engine
    if render_engine == "auto":
        for cand in ("BLENDER_EEVEE_NEXT", "BLENDER_EEVEE"):
            try:
                scene.render.engine = cand
                engine_used = cand
                break
            except TypeError:
                continue
        if engine_used == "auto":
            scene.render.engine = "CYCLES"
            engine_used = "CYCLES"
    else:
        scene.render.engine = render_engine
    if "EEVEE" in engine_used:
        try:
            scene.eevee.taa_render_samples = max(samples, 16)
        except AttributeError:
            pass
    else:
        scene.cycles.samples = samples
        scene.cycles.use_adaptive_sampling = True
        prefs = bpy.context.preferences.addons.get("cycles")
        if prefs:
            cprefs = prefs.preferences
            cprefs.compute_device_type = "NONE"
        scene.cycles.device = "CPU"
        try:
            scene.cycles.use_denoising = True
            scene.cycles.denoiser = "OPENIMAGEDENOISE"
        except (AttributeError, TypeError):
            pass

    scene.render.resolution_x = res[0]
    scene.render.resolution_y = res[1]
    scene.render.film_transparent = False
    try:
        scene.view_settings.view_transform = "Standard"   # true flat colors
    except TypeError:
        pass
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGB"
    return scene, engine_used


def add_sun(energy=0.35, direction=(0.4, -0.8, 1.0)):
    """Gentle sun for shading; world gives the flat white fill."""
    light = bpy.data.lights.new("sun", "SUN")
    light.energy = energy
    ob = bpy.data.objects.new("sun", light)
    bpy.context.scene.collection.objects.link(ob)
    d = Vector(direction).normalized()
    ob.rotation_euler = d.to_track_quat("Z", "Y").to_euler()
    return ob


# ---------------------------------------------------------------- meshes
def _read_stl(path):
    """Minimal binary/ASCII STL reader -> (verts Nx3 float, tris Nx3 int)."""
    with open(path, "rb") as f:
        blob = f.read()
    if blob[:5] == b"solid" and b"facet" in blob[:200]:
        verts, tris, cur = [], [], []
        for line in blob.decode("utf-8", "ignore").splitlines():
            p = line.split()
            if not p:
                continue
            if p[0] == "vertex":
                cur.append([float(v) for v in p[1:4]])
                if len(cur) == 3:
                    base = len(verts)
                    verts.extend(cur)
                    tris.append((base, base + 1, base + 2))
                    cur = []
        return np.array(verts), np.array(tris)
    n = struct.unpack("<I", blob[80:84])[0]
    rec = np.frombuffer(blob[84:84 + n * 50], dtype=np.uint8).reshape(n, 50)
    tri = rec[:, 12:48].copy().view("<f4").reshape(n, 3, 3).astype(np.float64)
    verts, tris, lut = [], [], {}
    for t in tri:
        idx = []
        for v in t:
            key = v.tobytes()
            if key not in lut:
                lut[key] = len(verts)
                verts.append(v)
            idx.append(lut[key])
        tris.append(idx)
    return np.array(verts), np.array(tris)


def import_stl(path, name, offset=None, location=(0, 0, 0), material=None):
    verts, tris = _read_stl(path)
    if offset is not None:
        verts = verts + np.asarray(offset)
    return _mesh_from_stl(verts, tris, name, location, material)


def import_stl_transformed(path, name, transform, material=None):
    """Import an STL and map vertices through `transform` (callable,
    (N,3) -> (N,3), in output units)."""
    verts, tris = _read_stl(path)
    verts = np.asarray(transform(np.asarray(verts, dtype=np.float64)))
    return _mesh_from_stl(verts, tris, name, (0, 0, 0), material)


def _mesh_from_stl(verts, tris, name, location, material):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts.tolist(), [], tris.tolist())
    mesh.update()
    ob = bpy.data.objects.new(name, mesh)
    bpy.context.scene.collection.objects.link(ob)
    ob.location = location
    if material:
        set_material(ob, material)
    return ob


def box(name, size, location, material, rounded=False):
    mesh = bpy.data.meshes.new(name)
    sx, sy, sz = size
    mesh.from_pydata(
        [(x * sx / 2, y * sy / 2, z * sz / 2) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)],
        [],
        [(0, 1, 3, 2), (4, 6, 7, 5), (0, 4, 5, 1), (2, 3, 7, 6), (0, 2, 6, 4), (1, 5, 7, 3)],
    )
    mesh.update()
    ob = bpy.data.objects.new(name, mesh)
    bpy.context.scene.collection.objects.link(ob)
    ob.location = location
    set_material(ob, material)
    return ob


def cylinder(name, radius, depth, location, material, axis="Z", verts=48):
    bpy.ops.mesh.primitive_cylinder_add(
        vertices=verts, radius=radius, depth=depth, location=location)
    ob = bpy.context.active_object
    ob.name = name
    if axis == "X":
        ob.rotation_euler = (0, math.pi / 2, 0)
    elif axis == "Y":
        ob.rotation_euler = (math.pi / 2, 0, 0)
    set_material(ob, material)
    return ob


# ---------------------------------------------------------------- curves
def poly_curve(name, points, radius=0.004, material=None, closed=False):
    """Straight polyline tube (curve with bevel). points: list of xyz."""
    curve = bpy.data.curves.new(name, "CURVE")
    curve.dimensions = "3D"
    spline = curve.splines.new("POLY")
    spline.points.add(len(points) - 1)
    for i, p in enumerate(points):
        spline.points[i].co = (p[0], p[1], p[2], 1.0)
    spline.use_cyclic_u = closed
    curve.bevel_depth = radius
    curve.bevel_resolution = 3
    curve.fill_mode = "FULL"
    ob = bpy.data.objects.new(name, curve)
    bpy.context.scene.collection.objects.link(ob)
    if material:
        set_material(ob, material)
    return ob


def set_polyline(ob, points):
    sp = ob.data.splines[0]
    assert len(sp.points) == len(points), f"{ob.name}: {len(sp.points)} vs {len(points)}"
    for i, p in enumerate(points):
        sp.points[i].co = (p[0], p[1], p[2], 1.0)


def keyframe_polyline(ob, frames, polylines, ease_boundary=True):
    """Keyframe every spline point at every frame (LINEAR by default prefs)."""
    n = len(frames)
    assert len(polylines) == n
    for k, pts in zip(frames, polylines):
        set_polyline(ob, pts)
        for pt in ob.data.splines[0].points:
            pt.keyframe_insert("co", frame=k)


def keyframe_linear(obj, data_path, frames, values, index=-1):
    for k, v in zip(frames, values):
        if index >= 0:
            getattr(obj, data_path)[index] = v
            obj.keyframe_insert(data_path, frame=k, index=index)
        else:
            setattr(obj, data_path, v)
            obj.keyframe_insert(data_path, frame=k)


def keyframe_mat_alpha(mat, frames, values):
    inp = mat.node_tree.nodes["Principled BSDF"].inputs["Alpha"]
    try:
        mat.blend_method = "BLEND"
    except (AttributeError, TypeError):
        pass
    for k, v in zip(frames, values):
        inp.default_value = v
        inp.keyframe_insert("default_value", frame=k)


# ---------------------------------------------------------------- materials
def make_material(name, color, alpha=1.0, emission=0.0, rough=1.0):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Alpha"].default_value = alpha
    try:
        bsdf.inputs["Roughness"].default_value = rough
        bsdf.inputs["Metallic"].default_value = 0.0
        bsdf.inputs["Specular IOR Level"].default_value = 0.05
    except KeyError:
        pass
    if emission > 0:
        bsdf.inputs["Emission Color"].default_value = (*color, 1.0)
        bsdf.inputs["Emission Strength"].default_value = emission
    if alpha < 1.0:
        try:
            mat.surface_render_method = "DITHERED"
        except AttributeError:
            mat.blend_method = "BLEND"
    return mat


def set_material(obj, mat):
    if isinstance(mat, (tuple, list)):
        mat = make_material(obj.name + "_mat", mat)
    if obj.data and hasattr(obj.data, "materials"):
        obj.data.materials.clear()
        obj.data.materials.append(mat)
    return mat


# ---------------------------------------------------------------- text
def load_font():
    for cand in (r"C:\Windows\Fonts\arialbd.ttf", r"C:\Windows\Fonts\arial.ttf"):
        try:
            f = bpy.data.fonts.load(cand)
            return f
        except RuntimeError:
            continue
    return None


def text_label(body, name, size, color=INK, align="LEFT"):
    curve = bpy.data.curves.new(name + "_font", "FONT")
    curve.body = body
    curve.size = size
    curve.align_x = align
    curve.align_y = "CENTER"
    f = load_font()
    if f:
        curve.font = f
    curve.extrude = size * 0.01
    ob = bpy.data.objects.new(name, curve)
    bpy.context.scene.collection.objects.link(ob)
    mat = make_material(name + "_ink", color)
    # Unlit text retains its intended contrast under the white studio lights.
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (0, 0, 0, 1)
    bsdf.inputs["Emission Color"].default_value = (*color, 1)
    bsdf.inputs["Emission Strength"].default_value = 1.0
    set_material(ob, mat)
    return ob


def track_camera(ob, camera):
    """Point a TEXT label's readable side (+Z) at the camera."""
    c = ob.constraints.new("COPY_ROTATION")
    c.target = camera
    c.owner_space = 'WORLD'
    c.target_space = 'WORLD'
    return c


# ---------------------------------------------------------------- arrows
def arrow(name, tail, head, radius, material, head_len_frac=0.28):
    """Cylinder shaft + cone head, grouped under an empty for easy animation."""
    tail = Vector(tail)
    head = Vector(head)
    d = head - tail
    L = d.length
    head_len = min(L * head_len_frac, 6 * radius)
    shaft_end = head - d.normalized() * head_len
    shaft = poly_curve(name + "_shaft", [tuple(tail), tuple(shaft_end)], radius=radius)
    set_material(shaft, material)
    cone = bpy.data.meshes.new(name + "_head")
    bpy.ops.mesh.primitive_cone_add(radius1=radius * 3.2, depth=head_len,
                                    vertices=24, location=(0, 0, 0))
    cone_ob = bpy.context.active_object
    cone_ob.name = name + "_head"
    cone_mesh = cone_ob.data
    bpy.data.objects.remove(cone_ob)
    cone_mesh.name = name + "_headmesh"
    head_ob = bpy.data.objects.new(name + "_head", cone_mesh)
    bpy.context.scene.collection.objects.link(head_ob)
    head_ob.location = head - d.normalized() * head_len / 2
    head_ob.rotation_euler = d.to_track_quat("Z", "Y").to_euler()
    set_material(head_ob, material)
    return shaft, head_ob


def aim_arrow(shaft, head_ob, tail, head):
    """Re-aim an existing arrow pair (from arrow()) at runtime."""
    tail = Vector(tail)
    head = Vector(head)
    d = head - tail
    L = d.length
    hl = min(L * 0.28, shaft.data.bevel_depth * 6)
    end = head - d.normalized() * hl
    set_polyline(shaft, [tuple(tail), tuple(end)])
    head_ob.location = head - d.normalized() * hl / 2
    head_ob.rotation_euler = d.to_track_quat("Z", "Y").to_euler()


# ---------------------------------------------------------------- camera
def look_rotation(location, target, up=(0, 1, 0)):
    """Explicit look-at rotation (camera looks along -Z, screen-up follows
    `up`). mathutils' to_track_quat assumes a Z-up world, which these
    anatomical-frame scenes are not, so the math is done by hand."""
    fwd = (Vector(target) - Vector(location)).normalized()
    up = Vector(up)
    right = fwd.cross(up)
    if right.length < 1e-9:                      # looking along `up`
        right = fwd.cross(Vector((0.0, 0.0, 1.0)))
    right.normalize()
    up2 = right.cross(fwd).normalized()
    # mathutils Matrix takes ROWS; camera basis needs these as COLUMNS
    return Matrix((right, up2, -fwd)).transposed().to_quaternion().to_euler()


def look_at(location, target, up=(0, 1, 0), lens=60):
    cam = bpy.data.cameras.new("cam")
    cam.lens = lens
    ob = bpy.data.objects.new("camera", cam)
    bpy.context.scene.collection.objects.link(ob)
    ob.rotation_euler = look_rotation(location, target, up)
    ob.location = location
    bpy.context.scene.camera = ob
    return ob


def render_frame(scene, frame, filepath):
    scene.frame_set(frame)
    scene.render.filepath = filepath
    bpy.ops.render.render(write_still=True)


def save_blend(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    bpy.ops.wm.save_as_mainfile(filepath=str(path), compress=True)
