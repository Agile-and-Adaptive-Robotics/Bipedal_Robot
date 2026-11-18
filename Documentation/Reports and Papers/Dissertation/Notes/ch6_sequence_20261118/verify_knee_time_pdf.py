"""Read every plotted vertex back from the PDF and compare with the saved MAT.

Run in the document runtime: python verify_knee_time_pdf.py FIGURE PLACED_PROOF OUTPUT_JSON
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
from pypdf import PdfReader
from pypdf.generic import ArrayObject


NUM = rb"[-+]?(?:\d*\.\d+|\d+)"
VERTEX = re.compile(rb"(" + NUM + rb") (" + NUM + rb") ([ml])(?=\s)")
WINDOW = re.compile(rb"n (" + NUM + rb") (" + NUM + rb") (" + NUM + rb") (" + NUM + rb") re W\* n")
FONT = re.compile(rb"/\S+ (" + NUM + rb") Tf")
T_RANGE = {"Main": (-1, 10), "Detail": (-0.020, 0.200)}
Y_RANGE = {"B": (-0.12, 2.34), "C": (-42.75, -41.40), "D": (-49.2, -37.8)}


def file_sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_curve(raw, key, t, y):
    group, row, signal = key.split("_")
    start = b"% TraceStart " + key.encode() + b"\n"
    end = b"% TraceEnd " + key.encode()
    previous, body = raw.split(start, 1)
    clip = WINDOW.findall(previous[-260:])[-1]
    l, b, w, h = [float(value) for value in clip]
    body = body.split(end, 1)[0]
    points = VERTEX.findall(body)
    assert points[0][2] == b"m" and all(p[2] == b"l" for p in points[1:])
    actual = np.array([(float(x), float(v)) for x, v, op in points])
    lo, hi = T_RANGE[group]
    ymin, ymax = Y_RANGE[row]
    if group == "Detail":
        valid = (t >= lo - 1e-10) & (t <= hi + 1e-10)
        t, y = t[valid], y[valid]
    model_count = len(t)
    if row == "B":
        # A rectangular pulse has two additional vertices for its vertical edges.
        jumps = np.flatnonzero(np.diff(y) != 0) + 1
        t = np.insert(t, jumps, t[jumps])
        y = np.insert(y, jumps, y[jumps - 1])
    expected = np.column_stack((l + w * (t - lo) / (hi - lo), b + h * (y - ymin) / (ymax - ymin)))
    assert actual.shape == expected.shape, (key, actual.shape, expected.shape)
    dx, dy = np.max(np.abs(actual - expected), axis=0)
    assert max(dx, dy) < 0.001, (key, dx, dy)
    zero = l - w * lo / (hi - lo)
    return {
        "source_samples": model_count,
        "recovered_pdf_vertices": len(points),
        "max_x_error_pt": dx,
        "max_y_error_pt": dy,
        "onset_x_pt": zero,
        "scale_pt_per_s": w / (hi - lo),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("figure", "proof", "output"):
        p.add_argument(arg, type=Path)
    p.add_argument(
        "--scientific-python", type=Path, default=Path("D:/Anaconda/envs/myo/python.exe")
    )
    p.add_argument("--proof-page", type=int, default=1)
    args = p.parse_args()
    page = PdfReader(args.figure).pages[0]
    # Compare the original decoded content streams with pdfTeX's embedded form.
    # Some pypdf releases add newlines when flattening an array in get_contents().
    # This figure already ends each source stream with its own delimiter.
    contents = page["/Contents"]
    if isinstance(contents, ArrayObject):
        raw = b"".join(ref.get_object().get_data() for ref in contents)
    else:
        raw = contents.get_object().get_data()
    with tempfile.TemporaryDirectory() as d:
        data = Path(d) / "arrays.npz"
        subprocess.run(
            [
                str(args.scientific_python),
                str(Path(__file__).with_name("make_knee_time_figure.py")),
                str(data),
                "--export-mat",
            ],
            check=True,
        )
        with np.load(data) as z:
            t = z["t"] - 30
            pulse = z["P"]
            ctrl = z["C"]
            curves = {}
            for name, signal in [
                ("B_stimulus", pulse[:, 8]),
                ("C_pulse", -np.rad2deg(pulse[:, 0])),
                ("C_control", -np.rad2deg(ctrl[:, 0])),
                ("D_extensor", pulse[:, 2]),
                ("D_flexor", pulse[:, 3]),
            ]:
                for group in T_RANGE:
                    key = group + "_" + name
                    curves[key] = checked_curve(raw, key, t, signal)
    for group in T_RANGE:
        axis_records = [v for k, v in curves.items() if k.startswith(group)]
        assert len({(v["onset_x_pt"], v["scale_pt_per_s"]) for v in axis_records}) == 1
    proof = PdfReader(args.proof).pages[args.proof_page - 1]
    proof_content = proof.get_contents().get_data()
    target = [
        k
        for k, v in proof["/Resources"]["/XObject"].items()
        if list(v.get_object().get("/BBox", [])) == list(page.mediabox)
    ]
    assert len(target) == 1, target
    embedded = proof["/Resources"]["/XObject"][target[0]]
    assert embedded.get_data() == raw
    schematic = page["/Resources"]["/XObject"]["/OriginalCircuit"].get_data()
    assert embedded["/Resources"]["/XObject"]["/OriginalCircuit"].get_data() == schematic
    before = proof_content.split(str(target[0]).encode() + b" Do")[0]
    transforms = re.findall(rb"(" + NUM + rb") 0 0 (" + NUM + rb") 0 0 cm", before[-170:])
    scale = float(transforms[-2][0])
    assert 0.9999 < scale < 1.0002
    fonts = [float(s) for s in FONT.findall(raw)]
    minimum = min(float(s) for s in FONT.findall(schematic)) * scale
    assert minimum >= 10 and min(fonts) * scale >= 10.4
    result = {
        "figure_sha256": file_sha(args.figure),
        "proof_sha256": file_sha(args.proof),
        "proof_page_one_based": args.proof_page,
        "matched_embedded_page_stream_sha256": hashlib.sha256(raw).hexdigest(),
        "matched_embedded_schematic_sha256": hashlib.sha256(schematic).hexdigest(),
        "latex_width_cm": 14,
        "placed_scale": scale,
        "minimum_retained_schematic_text_pt": minimum,
        "minimum_new_text_pt": min(fonts) * scale,
        "matched_vertices": curves,
        "max_pdf_coordinate_difference_pt": max(
            max(v["max_x_error_pt"], v["max_y_error_pt"]) for v in curves.values()
        ),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
