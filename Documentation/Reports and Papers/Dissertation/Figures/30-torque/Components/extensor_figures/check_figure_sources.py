"""Verify archived marker/curve positions, plotted data, fonts, and print rendering."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
from xml.etree import ElementTree

import numpy as np
import pymupdf
from PIL import Image
from scipy.io import loadmat


HERE = Path(__file__).resolve().parent
CHAPTER = HERE.parents[1]
REPO = next(p for p in HERE.parents if (p / "Testing_Data/2022_02_Festo").is_dir())
PREVIOUS = CHAPTER / "Deprecated/pinned_extensor_and_cam"
COLORS = {"hybrid": (1, .88, .25), "measured": (.25, .25, 1),
          "unoptimized": (.98, .53, .46), "optimized": (.8, .2, .71)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_previous(tests):
    result = []
    with pymupdf.open(PREVIOUS / "ExtPin_complete_publication.pdf") as document:
        draws = document[0].get_drawings()
    for j, trial in enumerate(tests):
        # Axis endpoints in the archived vector export, with ticks at -120 and +30
        # and the extra 15 degrees at the right that retained the final observations.
        left, right = ((100.75, 287.5), (348.25, 534.25), (595, 781.75))[j % 3]
        bottom = (258.75, 535.5)[j // 3]
        top = bottom - 195
        series = {}
        for name, target_color in COLORS.items():
            is_point = name in ("hybrid", "measured")
            original = []
            for draw in draws:
                color = draw.get("fill" if is_point else "color") or ()
                if tuple(round(x, 2) for x in color) != target_color:
                    continue
                rect = draw["rect"]
                if is_point:
                    center = (rect.tl + rect.br) * .5
                    if left - 4 <= center.x <= right + 4 and top - 4 <= center.y <= bottom + 4:
                        original.append(list(center))
                elif (len(draw["items"]) > 45 and rect.x0 >= left - .2 and rect.x1 <= right + .2
                      and rect.y0 >= top - .2 and rect.y1 <= bottom + .2):
                    original = [list(draw["items"][0][1])] + [list(item[2]) for item in draw["items"]]
            assert original, (j, name)
            old = np.array(original)
            if is_point:
                angles, torque = trial[name + "_angle"], trial[name]
                source = np.column_stack([left + (angles + 120) / 165 * (right - left),
                                          bottom - torque / 25 * 195])
                distances = np.sqrt(((old[:, None, :] - source[None, :, :]) ** 2).sum(axis=2))
                discrepancy = float(distances.min(axis=1).max())
                assert discrepancy < .006
                series[name] = {"source_points": len(angles), "old_visible_centers": len(old),
                                "maximum_position_error_pt": discrepancy}
            else:
                angles = (old[:, 0] - left) / (right - left) * 165 - 120
                torque = (bottom - old[:, 1]) / 195 * 25
                valid_nodes = torque > .2
                valid_values = np.isfinite(trial[name])
                measured = np.interp(angles[valid_nodes], trial["angle"][valid_values], trial[name][valid_values])
                discrepancy = float(np.max(np.abs(torque[valid_nodes] - measured)))
                assert discrepancy < .002
                series[name] = {"source_angle_samples": len(trial["angle"]),
                                "old_vertices_above_0.2_Nm": int(valid_nodes.sum()),
                                "maximum_curve_error_Nm": discrepancy}
        result.append({"panel": chr(65 + j), "original_test_id": int(trial["id"]),
                       "rest_cm": 100 * float(trial["rest"]), "series": series})
    return result


def check_torque_artwork(tests):
    """Read the final PDF paths and match all plotted markers and model vertices."""
    with pymupdf.open(CHAPTER / "pinned_extensor_torque.pdf") as doc:
        draws = doc[0].get_drawings()
        zero_ticks = [pymupdf.Rect(w[:4]) for w in doc[0].get_text("words") if w[4] == "0"]
    fills = {"hybrid": (239, 191, 22), "measured": (37, 63, 174)}
    strokes = {"unoptimized": (187, 77, 49), "optimized": (177, 42, 149)}
    result = []
    for j, trial in enumerate(tests):
        row, col = divmod(j, 3)
        left, top = 40 + 148 * col, 39 + 150 * row
        frame = pymupdf.Rect(left, top, left + 123, top + 92)
        inside = [d for d in draws if frame.contains(d["rect"].tl) and frame.contains(d["rect"].br)]
        # Frame bounds are in native PDF points. Read the signed marker position
        # in both dimensions; a point exactly on a clipped edge cannot pass.
        def project(angle, torque):
            return np.column_stack((left + (angle + 135) * 123 / 185,
                                    top + (25 - torque) * 92 / 27))

        def rgb(draw, key):
            return tuple(round(255 * value) for value in (draw[key] or ()))

        zero = project(0, 0)[0]
        rules = [d["rect"] for d in draws if rgb(d, "color") == (173, 178, 180)
                 and d["rect"].x0 == left and abs(d["rect"].y0 - zero[1]) < .0001]
        assert len(rules) == 1 and rules[0].x1 == frame.x1
        assert any(0 < left-tick.x1 < 8 and abs((tick.y0+tick.y1)/2-zero[1]) < .2
                   for tick in zero_ticks), "Missing signed torque origin tick"
        assert any(0 < tick.y0-frame.y1 < 6 and abs((tick.x0+tick.x1)/2-zero[0]) < .2
                   for tick in zero_ticks), "Missing angle origin tick"

        series = {}
        for name, color in fills.items():
            paths = [d for d in inside if rgb(d, "fill") == color]
            centers = np.array([list((d["rect"].tl + d["rect"].br) * .5) for d in paths])
            expected = project(trial[name + "_angle"], trial[name])
            assert centers.shape == expected.shape, (j, name, len(centers), len(expected))
            error = float(np.max(np.linalg.norm(centers - expected, axis=1)))
            assert error < .00025, (j, name, error)
            center_clearance = float(np.min([centers[:, 0]-frame.x0, frame.x1-centers[:, 0],
                                            centers[:, 1]-frame.y0, frame.y1-centers[:, 1]]))
            ink_clearance = min(min(d["rect"].x0-frame.x0, frame.x1-d["rect"].x1,
                                    d["rect"].y0-frame.y0, frame.y1-d["rect"].y1) - d["width"]/2
                                for d in paths)
            assert center_clearance > 5.5 and ink_clearance > 3.1
            series[name] = {"printed_markers": len(paths), "maximum_center_error_pt": error,
                            "minimum_center_clearance_pt": center_clearance,
                            "minimum_ink_clearance_pt": ink_clearance}
        for name, color in strokes.items():
            paths = [d for d in inside if rgb(d, "color") == color]
            assert len(paths) == 1
            segments = paths[0]["items"]
            vertices = np.array([list(segments[0][1])] + [list(segment[2]) for segment in segments])
            valid = np.isfinite(trial[name])
            expected = project(trial["angle"][valid], trial[name][valid])
            assert vertices.shape == expected.shape
            error = float(np.max(np.linalg.norm(vertices - expected, axis=1)))
            assert error < .0001, (j, name, error)
            series[name] = {"printed_finite_model_samples": len(vertices), "maximum_vertex_error_pt": error}
        result.append({"panel": chr(65 + j), "series": series})
    return {"common_angle_limits_deg": [-135, 50], "common_torque_limits_Nm": [-2, 25], "panels": result}


def check_published(stem, scratch):
    graphic = CHAPTER / (stem + ".pdf")
    process = subprocess.run(["pdftoppm", "-r", "216", "-png", "-singlefile", str(graphic),
                              str(scratch / stem)], check=True, capture_output=True)
    assert not process.stderr.strip(), process.stderr
    fonts = subprocess.run(["pdffonts", str(graphic)], check=True, capture_output=True).stdout.decode()
    assert all(line.split()[-5] == "yes" for line in fonts.splitlines()[2:] if line.strip())
    with Image.open(scratch / (stem + ".png")) as bitmap:
        numbers = np.asarray(bitmap.convert("RGB"))
        raster_size = list(bitmap.size)
        rgb_hash = hashlib.sha256(numbers.tobytes()).hexdigest()
    with pymupdf.open(graphic) as document:
        page = document[0]
        spans = [s for b in page.get_text("dict", flags=0)["blocks"]
                 for line in b.get("lines", []) for s in line["spans"] if s["text"].strip()]
        text = page.get_text()
        labels = (["Hybrid", "Measured", "Unoptimized", "Optimized", "(A)", "(B)", "(C)", "(D)", "(E)", "(F)"]
                  if stem == "pinned_extensor_torque" else
                  ["(A)", "(B)", "(C)", "Femur cylinder", "Femur wall", "Condyle envelope", "Lower tibia",
                   "Upper tibia", "Bend radius", "Optimized route", "Cam contacts", "Origin p1",
                   "Insertion p9", "Reference bones"])
        assert all(item in text for item in labels)
        assert min(s["size"] for s in spans) >= 9.99
        assert all("Arial" in s["font"] for s in spans)
        size = [page.rect.width, page.rect.height]
        font_range = [min(s["size"] for s in spans), max(s["size"] for s in spans)]
    svg = HERE / (stem + ".svg")
    tags = ElementTree.parse(svg).findall(".//{http://www.w3.org/2000/svg}text")
    assert len(tags) > 20
    return {"pdf_sha256": sha(graphic), "alt_sha256": sha(CHAPTER / (stem + ".alt.txt")),
            "editable_svg_sha256": sha(svg), "separate_svg_text_elements": len(tags),
            "native_width_height_pt": size, "native_type_size_range_pt": font_range,
            "native_216_dpi_render_px": raster_size, "rendered_rgb_sha256": rgb_hash, "pdffonts": fonts.splitlines()}


if __name__ == "__main__":
    data = loadmat(HERE / "extensor_data.mat", simplify_cells=True)
    raw = loadmat(REPO / "Testing_Data/2022_02_Festo/ExtPinBPASet.mat", simplify_cells=True)["ke"]
    for trial in data["tests"]:
        source = raw[int(trial["id"]) - 1]
        for field, original in (("measured_angle", "Aexp"), ("measured", "Mexp"),
                                ("hybrid_angle", "A_h"), ("hybrid", "M_h"),
                                ("angle", "Ak"), ("unoptimized", "M")):
            assert np.array_equal(trial[field], source[original], equal_nan=True)
    cam_source = REPO / "Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat"
    cam = loadmat(cam_source, simplify_cells=True)
    ctx, route = cam["ctxCamBest"], cam["predBest"]["routeInfo"]
    for i, pose in zip((0, -1), data["poses"]):
        transform = ctx["T_Pam"][:, :, i] @ ctx["T_ICR_t1"][:, :, i]
        xyz = route["raw"][:, :, i].copy()
        xyz[5:] = xyz[5:] @ transform[:3, :3].T + transform[:3, 3]
        assert np.allclose(pose["point"], xyz, atol=1e-14, rtol=0)
        assert np.array_equal(pose["active"], route["active"][:, i])
    assert sha(HERE / "cam_hardware_photo.jpg") == sha(PREVIOUS / "cam_extensor_hardware.jpg")
    with tempfile.TemporaryDirectory() as directory:
        result = {
            "graphics": {s: check_published(s, Path(directory)) for s in ("pinned_extensor_torque", "cam_extensor_route")},
            "pinned_tests": check_previous(data["tests"]),
            "pinned_final_artwork": check_torque_artwork(data["tests"]),
            "stored_pin_terms": data["terms"].tolist(),
            "active_cam_point_numbers": [(np.flatnonzero(p["active"]) + 1).tolist() for p in data["poses"]],
            "cam_angles_deg": [float(p["angle"]) for p in data["poses"]],
            "primary_inputs": {q.relative_to(REPO).as_posix(): sha(q) for q in (
                REPO / "Testing_Data/2022_02_Festo/ExtPinBPASet.mat",
                REPO / "Testing_Data/2022_02_Festo/minimizeExt10mmX3_results.mat", cam_source,
                REPO / "Code/Matlab/Bone_Mesh_Plots/Open_Sim_Bone_Geometry/Femur_Mesh_Points.txt",
                REPO / "Code/Matlab/Bone_Mesh_Plots/Open_Sim_Bone_Geometry/Tibia_Mesh_Points.txt",
                HERE / "cam_hardware_photo.jpg")},
        }
    (HERE / "verification.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

