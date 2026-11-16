"""Plot the native OpenSim preview and find cable centerline mesh crossings."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
import opensim as osim


def mesh_data(path):
    root = ET.parse(path)
    xyz = np.fromstring(root.findtext(".//Points/DataArray"), sep=" ").reshape(-1, 3)
    connectivity = np.fromstring(
        root.findtext(".//Polys/DataArray[@Name='connectivity']"), sep=" ", dtype=int
    )
    ends = np.fromstring(
        root.findtext(".//Polys/DataArray[@Name='offsets']"), sep=" ", dtype=int
    )
    faces = []
    for polygon in np.split(connectivity, ends[:-1]):
        faces.extend([polygon[0], polygon[i], polygon[i+1]] for i in range(1, len(polygon)-1))
    return xyz, np.asarray(faces)


def crossings(a, b, triangles):
    """Möller–Trumbore intersections, returned as fractions of a→b."""
    v, w = triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    h = np.cross(b-a, w)
    det = np.einsum("ij,ij->i", v, h)
    good = np.abs(det) > 1e-13
    inv = np.divide(1, det, out=np.zeros_like(det), where=good)
    s = a-triangles[:, 0]
    u = inv*np.einsum("ij,ij->i", s, h)
    q = np.cross(s, v)
    z = inv*(q @ (b-a))
    t = inv*np.einsum("ij,ij->i", w, q)
    valid = good & (u >= -1e-8) & (z >= -1e-8) & (u+z <= 1+1e-8)
    return np.unique(np.round(t[valid & (t > 1e-7) & (t < 1-1e-7)], 8))


def main(folder):
    model_file = folder / "tibialis_posterior_preview.osim"
    record = json.loads((folder / "tibialis_posterior_route.json").read_text())
    model = osim.Model(str(model_file))
    state = model.initSystem()
    model.realizePosition(state)
    bodies = model.getBodySet()
    pelvis = bodies.get("pelvis")

    def pelvis_point(body, point):
        p = bodies.get(body).findStationLocationInAnotherFrame(state, osim.Vec3(*point), pelvis)
        return [p.get(i) for i in range(3)]

    route = np.array([pelvis_point(b, p) for b, p in zip(record["bodies"], record["points_m"])])
    route -= route[1].copy()  # translate origin to the shank stop; retain pelvis orientation
    origin = np.array(pelvis_point(record["bodies"][1], record["points_m"][1]))
    meshes, local_meshes, findings = {}, {}, []
    xml = ET.parse(model_file)
    for body in ("tibia_r", "talus_r", "calcn_r", "toes_r"):
        for mesh in xml.findall(f".//Body[@name='{body}']/attached_geometry/Mesh"):
            assert mesh.findtext("socket_frame") == ".."
            file = Path(mesh.findtext("mesh_file"))
            xyz, faces = mesh_data(file)
            xyz *= np.fromstring(mesh.findtext("scale_factors"), sep=" ")
            local_meshes[file.stem] = (body, xyz[faces])
            xyz = np.array([pelvis_point(body, p) for p in xyz]) - origin
            meshes[file.stem] = (xyz, faces)
            for i, (a, b) in enumerate(zip(route[:-1], route[1:]), 1):
                ts = crossings(a, b, xyz[faces])
                if len(ts):
                    findings.append({"bone_mesh":file.stem, "segment":i,
                                     "locations_mm":(1000*(a+ts[:,None]*(b-a))).tolist()})

    fig, axes = plt.subplots(1, 2, figsize=(9, 6), layout="constrained", sharey=True)
    for ax, horizontal, label in zip(axes, (0, 2), ("Anterior +x (mm)", "Right / lateral +z (mm)")):
        for name, (xyz, faces) in meshes.items():
            color = "#b4bec8" if name != "fibula" else "#7c8b9e"
            ax.add_collection(PolyCollection(1000*xyz[faces][:, :, [horizontal, 1]],
                                             facecolor=color, edgecolor="none", alpha=.25))
        ax.plot(1000*route[:, horizontal], 1000*route[:, 1], "-o", color="#1262b2",
                linewidth=2, markersize=4, label="Cable centerline")
        ax.plot(0, 0, "ks", markersize=5, label="Shank stop")
        for hit in findings:
            loc = np.asarray(hit["locations_mm"])
            ax.scatter(loc[:, horizontal], loc[:, 1], s=32, c="#c73525", marker="x", zorder=5)
        ax.autoscale()
        ax.set_aspect("equal")
        ax.set_xlabel(label)
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Superior +y (mm)")
    axes[0].set_title("Sagittal view")
    axes[1].set_title("Frontal view")
    axes[1].legend(loc="upper right", fontsize=8)
    fig.suptitle("Right tibialis posterior | OpenSim default pose | pelvis axes")
    fig.savefig(folder/"tibialis_posterior_anatomical.png", dpi=200)
    fibula_offset = np.mean(meshes["fibula"][0], axis=0) - np.mean(meshes["tibia_r"][0], axis=0)
    rom_file = Path(__file__).resolve().parents[2] / "Routing27_v2" / "actuator_map_27.json"
    rom = json.loads(rom_file.read_text(encoding="utf-8"))["rom_table_deg"]
    sweep = []
    for coordinate in ("subtalar_angle_r", "ankle_angle_r"):
        joint = model.getCoordinateSet().get(coordinate)
        original = joint.getValue(state)
        hit_counts = {}
        angles = np.linspace(*rom[coordinate], 31)
        for angle in angles:
            joint.setValue(state, np.deg2rad(angle), False)
            model.realizePosition(state)
            for name, (body, triangles) in local_meshes.items():
                destination = bodies.get(body)
                points = []
                for owner, xyz in zip(record["bodies"], record["points_m"]):
                    point = bodies.get(owner).findStationLocationInAnotherFrame(
                        state, osim.Vec3(*xyz), destination)
                    points.append([point.get(k) for k in range(3)])
                for i, (a, b) in enumerate(zip(points[:-1], points[1:]), 1):
                    if len(crossings(np.array(a), np.array(b), triangles)):
                        hit_counts.setdefault(f"{name}:segment_{i}", []).append(float(angle))
        joint.setValue(state, original, False)
        sweep.append({"coordinate":coordinate, "angles_deg":angles.tolist(),
                      "mesh_crossings_deg":hit_counts})
    result = {"frame":"pelvis axes, translated to tibial stop", "pose":"OpenSim default",
              "fibula_minus_tibia_centroid_mm":(1000*fibula_offset).tolist(),
              "segments":["origin to shank stop", "shank stop to first foot guide",
                          "first foot guide to insertion"], "crossings":findings,
              "independent_coordinate_sweeps":sweep}
    (folder/"tibialis_posterior_mesh_check.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result))


if __name__ == "__main__":
    for directory in sys.argv[1:]:
        main(Path(directory))
