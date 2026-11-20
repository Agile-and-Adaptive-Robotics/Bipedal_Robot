"""Compare all plotted coordinates to their original numerical sources."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pymupdf
from scipy.interpolate import PchipInterpolator
from scipy.io import loadmat


HERE = Path(__file__).resolve().parent
DISS = HERE.parents[2]
REPO = DISS.parents[2]
ART = DISS / "Figures/30-torque"
OLD = ART / "Deprecated/physical_validation_20261119"
NEW_PANEL = (46.5, 5.625, 418.5, 298.875)
OLD_PANEL = (42.45, 2.10, 466.95, 239.10)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_pdf(path, *, panel, xlim, ylim, curves, points):
    left, top, right, bottom = panel
    with pymupdf.open(path) as doc:
        pg = doc[0]
        drawings = pg.get_drawings()
        page_size = list(pg.rect)[2:]
        type_sizes = sorted({round(s["size"], 3) for b in pg.get_text("dict")["blocks"]
                             for line in b.get("lines", []) for s in line["spans"]})
    checked = {}

    def compare(role, locations, data):
        xy = np.asarray(locations)
        x = xlim[0] + (xy[:, 0] - left) * (xlim[1] - xlim[0]) / (right - left)
        y = ylim[0] + (bottom - xy[:, 1]) * (ylim[1] - ylim[0]) / (bottom - top)
        ex = float(np.max(np.abs(x - data["x"])))
        ey = float(np.max(np.abs(y - data["y"])))
        assert ex < 0.001 and ey < 0.001, (role, ex, ey)
        checked[role] = dict(samples=len(x), max_degrees=ex, max_Nm=ey)

    for role, data in curves.items():
        paths = [d for d in drawings if len(d["items"]) == len(data["x"]) - 1
                 and (d["color"] == (0, 0, 0)) == (role == "human")]
        assert len(paths) == 1, role
        items = paths[0]["items"]
        compare(role, [items[0][1]] + [segment[2] for segment in items], data)
    counts = {}
    for role, data in points.items():
        scatters = [d for d in drawings if d["type"] == "f" and d["fill"]
                    and (d["fill"][2] > .98 and d["fill"][0] < .01 if role == "measured"
                         else d["fill"][0] > .98 and .8 < d["fill"][1] < .9)]
        count = len(data["x"])
        assert len(scatters) == count + 1  # the final marker belongs to the legend
        markers = scatters[:count]
        counts[role] = [s["seqno"] for s in markers]
        compare(role, [((s["rect"].x0 + s["rect"].x1) / 2,
                        (s["rect"].y0 + s["rect"].y1) / 2) for s in markers], data)
        expected_size = np.sqrt(data["area"])
        size_error = float(np.max(np.abs([s["rect"].width for s in markers] - expected_size)))
        assert size_error < .005, (role, size_error)
        opacity = [s["fill_opacity"] for s in markers]
        assert max(abs(a - data["alpha"]) for a in opacity) < .002
        checked[role].update(max_diameter_error_pt=size_error, opacity=opacity[0])
    assert min(counts["measured"]) > max(counts["prediction"])
    return dict(page_size_pt=page_size, panel_pt=panel, type_sizes_pt=type_sizes, series=checked)


def check_placed(proof):
    placed = {}
    with pymupdf.open(proof) as document:
        for page in document:
            text = page.get_text()
            for number, label in (("4.13", "Cam design of record, 620 kPa"),
                                  ("4.14", "Two-BPA design, 620 kPa")):
                if label not in text or f"Figure {number}:" not in text:
                    continue
                assert number not in placed
                spans = [s for b in page.get_text("dict")["blocks"] for line in b.get("lines", [])
                         for s in line["spans"]]
                legend_size = next(s["size"] for s in spans if s["text"] == label)
                curve = next(d for d in page.get_drawings() if len(d["items"]) == 99
                             and d["color"] and d["color"][0] > .75 and d["color"][2] > .65)
                scale = curve["width"] / 2.5
                assert abs(legend_size / 9.75 - scale) < .00001
                assert legend_size >= 10
                placed[number] = dict(pdf_page=page.number + 1, printed_page=text.splitlines()[0],
                                      legend_pt=legend_size, layout_scale=scale,
                                      plot_size_pt=[372 * scale, 293.25 * scale],
                                      caption=" ".join(text.split(f"Figure {number}:")[-1].split()))
        assert set(placed) == {"4.13", "4.14"}
        return dict(file=proof.resolve().relative_to(REPO).as_posix(), sha256=sha(proof),
                    total_pages=len(document), placed=placed)


def verify(active, output, proof):
    original_ext = OLD / "Ext20mm_physical_publication.pdf"
    original_flex = OLD / "Redesigned_flexor_pressure_publication.pdf"
    ext = ART / "extensor_validation.pdf"
    assert sha(ext) == sha(original_ext)
    source = json.loads((HERE / "flexor_data.json").read_text(encoding="utf-8"))
    series = {s["name"]: s for s in source["series"]}
    curves = {k: series[k] for k in ("design", "human")}
    points = {k: series[k] for k in ("measured", "prediction")}
    kwargs = dict(ylim=source["ylim"], curves=curves, points=points)
    old_flex = compare_pdf(original_flex, panel=OLD_PANEL, xlim=source["xlim"], **kwargs)
    new_flex = compare_pdf(active, panel=NEW_PANEL, xlim=(-125, 15), **kwargs)
    assert new_flex["page_size_pt"] == [422, 338]
    reference = REPO / "Testing_Data/2022_02_Festo/OpenSim_Bifem_Results.txt"
    human = np.loadtxt(reference, skiprows=7)
    assert np.array_equal(human[:, 3], series["human"]["y"])
    assert np.array_equal(human[:, 1], series["human"]["x"])
    trial = points["measured"]
    margin = np.abs(trial["y"]) - np.abs(PchipInterpolator(human[:, 1], human[:, 3])(trial["x"]))
    model = REPO / "Code/Matlab/Mesh_Optimization/Results/Bifemsh_20mm_Result.mat"
    rest_mm = float(loadmat(model, variable_names=["xBest"], squeeze_me=True)["xBest"][6])*1000
    flex_pred = ART / "Components/publication_20261028/flexor_pressure_publication.json"
    evaluated = json.loads(flex_pred.read_text(encoding="utf-8"))
    assert np.array_equal(evaluated["prediction"], points["prediction"]["y"])
    assert np.array_equal(evaluated["angle"], trial["x"])
    ext_csv = HERE / "extensor_trials.csv"
    with ext_csv.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    def column(key):
        return np.array([float(row[key]) for row in rows])
    ext_points = {role: dict(x=column("K_deg"), y=column(key),
                             area=60*column("Pressure_kPa")/620,
                             alpha=.75 if role == "measured" else 1)
                  for role, key in (("measured", "measured_torque_Nm_RECORD"),
                                    ("prediction", "pred_at_measured_P_Nm"))}
    cam = model.with_name("Vas_Pam_20mm_Result_cam_20261014_0416.mat")
    design = loadmat(cam, variable_names=["predBest", "ctxCamBest"], simplify_cells=True)
    ctx = design["ctxCamBest"]
    ext_curves = {"design": dict(x=ctx["phiD"], y=design["predBest"]["TorqueZ"]),
                  "human": dict(x=ctx["humanAngleD"], y=ctx["humanTorque"])}
    ext_check = compare_pdf(ext, panel=NEW_PANEL, xlim=(-125, 15), ylim=(20, 110),
                            curves=ext_curves, points=ext_points)
    old_fig_sha = "2f5bc4ea2e6df965338db82762b148660c31d9afc949a0d4bad3f0c9cdb5fb96"
    moved = ART / "Components/flexor_comparison.fig"
    assert sha(moved) == old_fig_sha
    source_figure = ART / "Components/figure_polish_20261030/Redesigned_flexor_pressure_publication.fig"
    inputs = [original_ext, original_flex, OLD / "Redesigned_flexor_pressure_publication.alt.txt",
              ext, active, ART / "extensor_validation.alt.txt", ART / "flexor_validation.alt.txt",
              source_figure, model, cam, reference, flex_pred, ext_csv, moved,
              ART / "Components/publication_20261028/publication_flexor_pressure.m",
              HERE / "read_flexor.m", HERE / "flexor_data.json", HERE / "draw_flexor.py"]
    result = dict(
        old_flexor=old_flex, final_flexor=new_flex, unchanged_extensor=ext_check,
        flexor_margin_Nm=dict(low=float(margin.min()), low_degrees=trial["x"][np.argmin(margin)],
                              high=float(margin.max()), high_degrees=trial["x"][np.argmax(margin)]),
        extensor_margin_Nm=[column("meas_minus_human").min(), column("meas_minus_human").max()],
        model_rest_mm=rest_mm,
        files={p.resolve().relative_to(REPO).as_posix():sha(p) for p in inputs})
    if proof is not None:
        result["combined_proof"] = check_placed(proof)
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flexor", type=Path, default=ART / "flexor_validation.pdf")
    parser.add_argument("--out", type=Path, default=HERE / "verification.json")
    parser.add_argument("--proof", type=Path, help="Verify both complete printed figure pages.")
    args = parser.parse_args()
    verify(args.flexor, args.out, args.proof)
