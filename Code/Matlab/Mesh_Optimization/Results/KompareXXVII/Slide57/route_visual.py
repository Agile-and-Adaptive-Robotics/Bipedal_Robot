"""Draw the replaceable right-hand image for defense slide 57."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if p.name == "Bipedal_Robot")
ROUTING = ROOT / "Code/Matlab/Mesh_Optimization/Functional_Routing_20261027"
GEOMETRY = ROUTING / "geometry_20261103/extraction_simbody.json"
MODEL = ROOT / "Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim"
SOLEUS = (HERE.parents[1] / "BiPulley_surveys"
          / "BiPulley_survey2_soleus_r_nB2_G1_20261023_052012.mat")
COLORS = {"single": "#157349", "bundle": "#1958AA", "tackle": "#BD520F"}


def project(points):
    p = np.atleast_2d(points)
    return np.c_[.80 * p[:, 0] - .60 * p[:, 2], p[:, 1]]


def line(ax, points, color, width=1.65, style="-", alpha=1.):
    p = project(points)
    ax.plot(*p.T, color=color, lw=width, linestyle=style, alpha=alpha)


def group(row):
    if max(int(row[key]) for key in ("tackleLineParts_c1", "tackleLineParts_c2")) > 1:
        return "tackle"
    return "bundle" if int(row["nBPA"]) > 1 else "single"


def relative_frame(neutral, body, parent):
    child, reference = neutral[body], neutral[parent]
    inverse_rotation = np.asarray(reference["R"]).T
    transform = np.eye(4)
    transform[:3, :3] = inverse_rotation @ child["R"]
    transform[:3, 3] = inverse_rotation @ (np.asarray(child["p"]) - reference["p"])
    return transform


def draw_skeleton(ax, survey, neutral, geometry_dir):
    sys.path.insert(0, str(ROUTING / "campaign_20261103"))
    from fr_plotPilotAnatomy import mesh_data

    foot = loadmat(SOLEUS, squeeze_me=True, struct_as_record=False)
    foot_grid = foot["T"]
    ankle_index, subtalar_index = [
        np.argmin(np.abs(np.linspace(start, end, foot_grid.shape[2])))
        for start, end in foot["spec"].kin.thetaRangesDeg
    ]
    ankle = foot_grid[:, :, ankle_index, 0]
    subtalar = foot_grid[:, :, subtalar_index, 1]
    np.testing.assert_allclose(ankle[:3, 3], survey["ankP"], atol=1e-9)
    np.testing.assert_allclose((ankle @ subtalar)[:3, 3], survey["subP"], atol=1e-9)
    talus = survey["Ttib"] @ ankle
    calcaneus = talus @ subtalar
    frames = {
        "torso": relative_frame(neutral, "torso", "pelvis"),
        "pelvis": np.eye(4),
        "femur_r": survey["Thip"],
        "tibia_r": survey["Ttib"],
        "talus_r": talus,
        "calcn_r": calcaneus,
        "toes_r": calcaneus @ relative_frame(neutral, "toes_r", "calcn_r"),
    }
    model = ET.parse(MODEL)
    surfaces = []
    for body, transform in frames.items():
        for mesh in model.findall(f".//Body[@name='{body}']/attached_geometry/Mesh"):
            vertices, faces = mesh_data(geometry_dir / mesh.findtext("mesh_file"))
            vertices *= np.fromstring(mesh.findtext("scale_factors"), sep=" ")
            vertices = vertices @ transform[:3, :3].T + transform[:3, 3]
            surfaces.append(vertices[faces])
    triangles = np.concatenate(surfaces)
    depth_axis = [.60, 0, .80]
    order = np.argsort(triangles.mean(axis=1) @ depth_axis)
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    norm = np.maximum(np.linalg.norm(normals, axis=1), np.finfo(float).eps)
    shade = np.abs(normals @ [.28, .76, .59]) / norm
    facecolors = np.clip(np.array([.72, .75, .77]) + .22 * shade[:, None], 0, 1)
    xy = project(triangles[order].reshape(-1, 3)).reshape(-1, 3, 2)
    ax.add_collection(PolyCollection(xy, facecolors=facecolors[order], edgecolors="none",
                                     antialiased=False, zorder=0))
    ax.autoscale_view()


def connected_psoas(raw, iliacus, torso_to_pelvis):
    proximal = np.asarray(raw[:2]) @ torso_to_pelvis[:3, :3].T + torso_to_pelvis[:3, 3]
    np.testing.assert_allclose(proximal[1], iliacus[1], atol=1e-9)
    # Draw the shared distal connection at the saved iliacus attachment.
    return np.vstack([proximal, iliacus[-1]])


def draw_erector(ax, record):
    points = np.array(record["points_m"], dtype=float)
    transform = np.array(record["torso_to_pelvis"])
    for i, frame in enumerate(record["frames"]):
        if frame == "torso":
            points[i] = transform[:3, :3] @ points[i] + transform[:3, 3]
        else:
            assert frame == "pelvis"
    assert points.shape == (2, 3)
    line(ax, points, COLORS["bundle"])
    ax.scatter(*project(points).T, s=7, color=COLORS["bundle"], zorder=3)


def draw_table(selected, record, stem):
    fig, ax = plt.subplots(figsize=(13.04, 5.63), facecolor="white")
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax.axis("off")
    labels = {
        "bifemlh_r": "Biceps Femoris, Long Head",
        "bifemsh_r": "Biceps Femoris, Short Head",
        "med_gas_r": "Medial Gastrocnemius",
        "lat_gas_r": "Lateral Gastrocnemius",
        "tib_ant_r": "Tibialis Anterior",
        "psoas_r": "Psoas",
    }
    for start, column in ((0, selected[:14]), (.5, selected[14:])):
        xs = (start + .025, start + .35, start + .44)
        for x, label, align in zip(xs, ("Actuator site", "BPAs", "Parts"),
                                   ("left", "center", "center"), strict=True):
            ax.text(x, .983, label, ha=align, va="top", fontsize=18, weight="bold")
        for i, row in enumerate(column):
            name = row["site"]
            count = int(row["nBPA"])
            parts = max(int(row[key]) for key in ("tackleLineParts_c1", "tackleLineParts_c2"))
            color = COLORS[group(row)]
            if name == "ercspn_r":
                count, parts = record["count"], record["gain"]
                color = COLORS["bundle"]
            y = .91 - i * .061
            ax.text(xs[0], y, labels.get(name, row["muscle"]), color=color,
                    va="center", fontsize=15.5)
            for x, value in zip(xs[1:], (count, parts), strict=True):
                ax.text(x, y, str(value), ha="center", va="center", fontsize=15.5)
    fig.savefig(HERE / f"{stem}.png", dpi=240)
    fig.savefig(HERE / f"{stem}.pdf")
    plt.close(fig)


def render(branch, geometry_dir, schematic_only):
    with (HERE.parent / "route_configurations.csv").open(encoding="utf-8-sig", newline="") as f:
        selected = [r for r in csv.DictReader(f) if r["side"] == "selected"]
    assert len(selected) == 27
    m = loadmat(HERE / "route_projection.mat", squeeze_me=True, struct_as_record=False)
    assert all(x.site == r["site"] for x, r in zip(m["sel"], selected, strict=True))
    rows = {r["site"]: np.asarray(x.P) for r, x in zip(selected, m["Pel"], strict=True)}
    neutral = json.loads(GEOMETRY.read_text(encoding="utf-8"))["bodyTransformsInGround"]["neutral"]
    psoas = connected_psoas(rows["psoas_r"], rows["iliacus_r"],
                            relative_frame(neutral, "torso", "pelvis"))
    assert np.array_equal(psoas[-1], rows["iliacus_r"][-1])
    record = json.loads((HERE / f"erector_rerun/erector_{branch}.json").read_text())

    plt.rcParams.update({"font.family": "Arial"})
    fig = plt.figure(figsize=(5.8, 5.5), facecolor="white")
    ax = fig.add_axes([.0, .02, .73, .96])
    draw_skeleton(ax, m, neutral, geometry_dir)
    for row in selected:
        name = row["site"]
        if name in ("psoas_r", "ercspn_r"):
            continue
        mode = group(row)
        line(ax, rows[name], COLORS[mode], style="--" if mode == "tackle" else "-")
        ax.scatter(*project(rows[name][[0, -1]]).T, s=7, color=COLORS[mode], zorder=3)
    line(ax, psoas[:2], COLORS["bundle"])
    ax.scatter(*project(psoas[[0]]).T, s=7, color=COLORS["bundle"], zorder=3)
    draw_erector(ax, record)
    ax.set(aspect="equal")
    ax.axis("off")
    suffix = "" if branch == "compliant" else "_rigid"
    fig.savefig(HERE / f"route_schematic{suffix}.png", dpi=300)
    fig.savefig(HERE / f"route_schematic{suffix}.pdf")
    plt.close(fig)
    if not schematic_only:
        draw_table(selected, record, f"actuator_site_table{suffix}")
    counts = [record["count"] if row["site"] == "ercspn_r" else int(row["nBPA"]) for row in selected]
    print(f"{sum(counts)} BPAs across {len(counts)} sites, counts {dict(sorted(Counter(counts).items()))}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch", choices=("rigid", "compliant"), default="compliant")
    parser.add_argument("--geometry", type=Path, default=Path("C:/OpenSim 4.6/Geometry"))
    parser.add_argument("--schematic-only", action="store_true")
    args = parser.parse_args()
    render(args.branch, args.geometry, args.schematic_only)
