"""Build Figure 6.4 from the matched 1 ms MATLAB replay. Run with document Python.

Example: python make_knee_time_figure.py OUTPUT_DIR [--replay path/to/reader_summary.mat]
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / "Code/Matlab/SNS_Simscape").is_dir())
DISS = ROOT / "Documentation/Reports and Papers/Dissertation"
MATDIR = ROOT / "Code/Matlab/SNS_Simscape/results/knee_transient_reflex_20261028"
OLD = DISS / "Figures/60-sns/Deprecated/20261118_knee_time_alignment"
STEM = "sns_knee_reflex_authoritative"
DATA_HASH = "01a67c749c0a317986cbcc7d55192f0795770c722b15bd4ee52f89cecfd9309c"
OLD_HASH = "771c02ec1777d181ee19d36d69b4dac635ffc2891e2b55be64adf836480a3747"
W, H = 14 * 72 / 2.54, 472.0
MAIN = (51.0, 233.0, -1.0, 10.0)
DETAIL = (263.0, 384.0, -0.020, 0.200)
ROWS = {
    "B": (237.0, 46.0, -0.12, 2.34),
    "C": (145.0, 65.0, -42.75, -41.40),
    "D": (53.0, 65.0, -49.2, -37.8),
}
BLUE, RED, GREY = "#4477AA", "#EE6677", "#484848"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def export_mat(filename, replay):
    # SciPy is supplied by the project scientific runtime.
    from scipy.io import loadmat

    raw = (MATDIR / "replay_summary.mat").read_bytes()
    assert sha(raw) == DATA_HASH
    v = loadmat(io.BytesIO(raw), simplify_cells=True)
    names = ["th", "thd", "V_MN_ext", "V_MN_flex", "A_ext", "A_flex", "F_ext", "F_flex", "stim"]
    assert list(v["names"]) == names
    if replay:
        check = loadmat(replay, simplify_cells=True)
        assert all(np.array_equal(v[k], check[k]) for k in ("t", "C", "P", "F", "names"))
    np.savez(filename, **{k: v[k] for k in ("t", "C", "P", "F")})


def xy(col, row, t, value):
    l, r, lo, hi = col
    b, h, ymin, ymax = ROWS[row]
    return l + (t - lo) * (r - l) / (hi - lo), b + (value - ymin) * h / (ymax - ymin)


def text(c, x, y, s, bold=False, size=10.5, align="left", color="#181818"):
    c.setFillColor(color)
    c.setFont("KneeArialBold" if bold else "KneeArial", size)
    {"left": c.drawString, "center": c.drawCentredString, "right": c.drawRightString}[align](
        x, y, s
    )


def stroke(c, x0, y0, x1, y1, color, width=0.6, dash=()):
    c.setLineWidth(width)
    c.setStrokeColor(color)
    c.setDash(dash)
    c.line(x0, y0, x1, y1)


def background(c, row, col, ticks, yticks):
    l, r, lo, hi = col
    b, h, ymin, ymax = ROWS[row]
    for tick in ticks:
        x, _ = xy(col, row, tick, ymin)
        stroke(c, x, b, x, b + h, "#dedede", 0.45)
        stroke(c, x, b, x, b - 2.6, "#444444", 0.65)
        if row == "D":
            display = tick * 1000 if col == DETAIL else tick
            text(c, x, b - 14, f"{display:.0f}", size=10.4, align="center")
    for tick in yticks:
        _, y = xy(col, row, lo, tick)
        stroke(c, l, y, r, y, "#dedede", 0.45)
        if col == MAIN:
            stroke(c, l - 2.8, y, l, y, "#444444", 0.65)
            text(c, l - 5, y - 3.4, f"{tick:g}".replace("-", "−"), size=10.4, align="right")
    stroke(c, l, b, l, b + h, "#4b4b4b", 0.7)
    stroke(c, l, b, r, b, "#4b4b4b", 0.7)
    onset, _ = xy(col, row, 0, ymin)
    stroke(c, onset, b, onset, b + h, "#4b4b4b", 0.75, (1.4, 1.5))


def series(c, row, col, t, values, label, color, width, dash=(), step=False):
    l, r, lo, hi = col
    b, h, ymin, ymax = ROWS[row]
    # All 50,001 raw measurements are present in each main vector path. PDF
    # clipping selects the settled response; the onset detail uses its raw subset.
    use = (
        np.arange(len(t)) if col == MAIN else np.flatnonzero((t >= lo - 1e-10) & (t <= hi + 1e-10))
    )
    tx, vy = xy(col, row, t[use], values[use])
    c.saveState()
    bounds = c.beginPath()
    bounds.rect(l, b, r - l, h)
    c.clipPath(bounds, stroke=0, fill=0)
    c.setStrokeColor(color)
    c.setLineWidth(width)
    c.setDash(dash)
    c.setLineJoin(1)
    c.setLineCap(0)
    key = ("Main" if col == MAIN else "Detail") + "_" + row + "_" + label
    c._code.append("% TraceStart " + key)
    path = c.beginPath()
    path.moveTo(tx[0], vy[0])
    last = vy[0]
    jumps = 0
    for x, y in zip(tx[1:], vy[1:]):
        if step and y != last:
            path.lineTo(x, last)
            jumps += 1
        path.lineTo(x, y)
        last = y
    c.drawPath(path, stroke=1, fill=0)
    c._code.append("% TraceEnd " + key)
    c.restoreState()
    return key, {
        "samples": len(use),
        "vertices": len(use) + jumps,
        "raw_source_span_s": [float(t[use[0]]), float(t[use[-1]])],
        "visible_samples": int(np.count_nonzero((t >= lo - 1e-10) & (t <= hi + 1e-10))),
    }


def sample_label(c, x, y, title, color, dash=()):
    stroke(c, x, y + 3, x + 12, y + 3, color, 1.35, dash)
    text(c, x + 15, y, title, size=10.4)


def render(z):
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.pdfgen import canvas

    for label, name in [("KneeArial", "arial.ttf"), ("KneeArialBold", "arialbd.ttf")]:
        pdfmetrics.registerFont(TTFont(label, str(Path("C:/Windows/Fonts") / name)))
    t = z["t"] - 30
    p = z["P"]
    ctrl = z["C"]
    angle = -np.rad2deg(p[:, 0])
    base_angle = -np.rad2deg(ctrl[:, 0])
    refs = {}
    out = io.BytesIO()
    c = canvas.Canvas(out, pagesize=(W, H), invariant=True, pageCompression=1)
    c.setTitle("Knee reflex response on aligned time axes")
    # A retains the exact vector drawing and all six connections from the source.
    c.setFillColor("#f4f5f6")
    c.roundRect(249, 7, 145, H - 163, 4, stroke=0, fill=1)
    text(c, (MAIN[0] + MAIN[1]) / 2, H - 164, "Response and recovery", True, 11, align="center")
    text(c, (DETAIL[0] + DETAIL[1]) / 2, H - 164, "Pulse-onset detail", True, 11, align="center")
    for row, y, name in [
        ("B", 292, "Injected extensor Ia current (nA)"),
        ("C", 221, "Knee angle (deg)"),
        ("D", 129, "MN voltage (mV)"),
    ]:
        text(c, 8, y, row, True, 11.5)
        text(c, 25, y, name, True, 10.6)
    sample_label(c, 233, 221, "Stimulus", BLUE)
    sample_label(c, 298, 221, "No stimulus", GREY, (4, 2.4))
    sample_label(c, 233, 129, "Extensor", BLUE)
    sample_label(c, 298, 129, "Flexor", RED, (4, 2.4))
    for col in (MAIN, DETAIL):
        ticks = [0, 2, 4, 6, 8, 10] if col == MAIN else [0, 0.05, 0.10, 0.15, 0.20]
        for row, vals in [("B", [0, 2]), ("C", [-42.5, -42, -41.5]), ("D", [-48, -44, -40])]:
            background(c, row, col, ticks, vals)
        for row, v, lab, color, w, dash, step in [
            ("B", p[:, 8], "stimulus", "#232323", 1.35, (), True),
            ("C", angle, "pulse", BLUE, 1.35, (), False),
            ("C", base_angle, "control", GREY, 1.05, (4, 2.4), False),
            ("D", p[:, 2], "extensor", BLUE, 1.25, (), False),
            ("D", p[:, 3], "flexor", RED, 1.2, (4, 2.4), False),
        ]:
            k, entry = series(c, row, col, t, v, lab, color, w, dash, step)
            refs[k] = entry
    text(c, (MAIN[0] + MAIN[1]) / 2, 18, "Time from pulse onset (s)", size=10.7, align="center")
    text(c, (DETAIL[0] + DETAIL[1]) / 2, 18, "Time from onset (ms)", size=10.7, align="center")
    # Exact annotated samples. The knee peak is relative to concurrent control.
    ipeak = int(np.flatnonzero(t >= 0)[np.argmax((angle - base_angle)[t >= 0])])
    px, py = xy(DETAIL, "C", t[ipeak], angle[ipeak])
    c.setStrokeColor(BLUE)
    c.setFillColor("white")
    c.setLineWidth(0.8)
    c.setDash([])
    c.circle(px, py, 1.8, stroke=1, fill=1)
    stroke(c, px + 1.4, py - 1.4, px + 8, py - 22, "#494949", 0.6)
    text(c, px + 9, py - 33, "Peak +1.09°", size=10.4, align="center")
    text(c, px + 9, py - 45, "at 111 ms", size=10.4, align="center")
    irec = int(np.argmin(abs(t - 6.355)))
    rx, ry = xy(MAIN, "C", t[irec], angle[irec])
    stroke(c, rx, ry + 1, rx + 3, ry + 16, "#494949", 0.6)
    text(c, rx + 1, ry + 31, "Recovery", size=10.4, align="center")
    text(c, rx + 1, ry + 20, "6.355 s", size=10.4, align="center")
    text(c, 330, 270, "2 nA for 30 ms", size=10.4, align="center")
    c.showPage()
    c.save()
    return out.getvalue(), refs


def add_circuit(plots, source):
    from pypdf import PdfReader, PdfWriter
    from pypdf.generic import ArrayObject, DecodedStreamObject, FloatObject, NameObject

    assert sha(source) == OLD_HASH
    old = PdfReader(io.BytesIO(source)).pages[0]
    raw = old.get_contents().get_data()
    boundary = b"/A1 gs 1 g 0 j 0 w 0 G 1 g\n\n63.496063 216.850394 m"
    assert raw.count(boundary) == 1
    drawing = raw[: raw.index(boundary)]
    writer = PdfWriter(clone_from=PdfReader(io.BytesIO(plots)))
    page = writer.pages[0]
    form = DecodedStreamObject()
    form.set_data(drawing)
    form.update(
        {
            NameObject("/Type"): NameObject("/XObject"),
            NameObject("/Subtype"): NameObject("/Form"),
            NameObject("/BBox"): ArrayObject([FloatObject(v) for v in (0, 274.5, W, 425.19685)]),
            NameObject("/Resources"): old["/Resources"].clone(writer),
        }
    )
    from pypdf.generic import DictionaryObject

    resources = page["/Resources"]
    xobjects = resources.setdefault(NameObject("/XObject"), DictionaryObject())
    xobjects[NameObject("/OriginalCircuit")] = writer._add_object(form.flate_encode())
    overlay = DecodedStreamObject()
    overlay.set_data(f"q 1 0 0 1 0 {H-425.19685:.8f} cm /OriginalCircuit Do Q\n".encode())
    page[NameObject("/Contents")] = ArrayObject(
        [page.raw_get("/Contents"), writer._add_object(overlay)]
    )
    output = io.BytesIO()
    writer.write(output)
    pdf = output.getvalue()
    verify = PdfReader(io.BytesIO(pdf)).pages[0]
    assert verify["/Resources"]["/XObject"]["/OriginalCircuit"].get_data() == drawing
    return pdf, sha(drawing), verify.get_contents().get_data()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--replay", type=Path)
    parser.add_argument(
        "--scientific-python", type=Path, default=Path("D:/Anaconda/envs/myo/python.exe")
    )
    parser.add_argument("--export-mat", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.export_mat:
        export_mat(args.output, args.replay)
        return
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        data = Path(tmp) / "arrays.npz"
        command = [
            str(args.scientific_python),
            str(Path(__file__).resolve()),
            str(data),
            "--export-mat",
        ]
        if args.replay:
            command.extend(["--replay", str(args.replay)])
        subprocess.run(command, check=True)
        with np.load(data) as z:
            plots, records = render(z)
            t, p, c = z["t"], z["P"], z["C"]
            relative = t - 30
            extrema = np.flatnonzero(p[:, 8] > 0)
            assert (
                t.size == 50001
                and p.shape == c.shape == (50001, 9)
                and np.allclose(np.diff(t), 0.001, atol=1e-13, rtol=0)
            )
            assert (
                np.array_equal(extrema, np.arange(30000, 30030))
                and np.all(p[extrema, 8] == 2)
                and np.all(c[:, 8] == 0)
            )
            baseline = np.mean(c[(t >= 25) & (t < 30), 0])
            extension = -np.rad2deg(p[:, 0] - c[:, 0])
            pick = 30000 + np.argmax(extension[30000:])
            bad = np.flatnonzero(
                (t >= 30.03)
                & ((abs(np.rad2deg(p[:, 0] - baseline)) > 0.01) | (abs(p[:, 1]) > 1e-4))
            )
            summary = json.loads((MATDIR / "replay_summary.json").read_text())
            for got, exp in [
                (relative[pick], summary["peak_latency_s"]),
                (relative[bad[-1] + 1], summary["recovery_s"]),
                (extension[pick], summary["paper_extension_excursion_deg"]),
            ]:
                assert np.isclose(got, exp, atol=1e-12, rtol=0)
            numeric = {
                "peak_latency_s": relative[pick],
                "peak_extension_deg": extension[pick],
                "recovery_s": relative[bad[-1] + 1],
                "control_baseline_deg": -np.rad2deg(baseline),
                "pulse_rows": int(len(extrema)),
                "array_sha256_C_order": {
                    k: sha(np.ascontiguousarray(z[k]).tobytes()) for k in ("t", "P", "C", "F")
                },
            }
    old = OLD / (STEM + ".pdf")
    if not old.exists():
        old = DISS / "Notes/ch6_authoritative_20261116/CAD_sns_compact" / (STEM + ".pdf")
    pdf, schematic_sha, content = add_circuit(plots, old.read_bytes())
    (args.output / (STEM + ".pdf")).write_bytes(pdf)
    for key, v in records.items():
        token = key.encode()
        part = content.split(b"% TraceStart " + token + b"\n", 1)[1].split(
            b"% TraceEnd " + token, 1
        )[0]
        assert len(re.findall(rb"(?<!\S)[ml](?!\S)", part)) == v["vertices"], key
    alt = (
        "A shows the six antagonist knee-reflex connections. B, C, and D share time zero at the start of a 2 nA, 30 ms extensor Ia pulse. "
        "Their left axes cover −1 to +10 s. Matching right-hand details enlarge −20 to +200 ms. B shows injected current. "
        "C shows knee angles in the stimulated and no-stimulus trials, with flexion negative. The knee extended 1.091° above the concurrent control "
        "at 111 ms and met the angle and speed recovery criteria at 6.355 s. D shows extensor and flexor motor-neuron membrane voltages in the stimulated trial.\n"
    )
    (args.output / (STEM + ".alt.txt")).write_text(alt, encoding="utf-8")
    report = {
        "source_mat": str((MATDIR / "replay_summary.mat").relative_to(ROOT)),
        "source_mat_sha256": DATA_HASH,
        "source_summary_sha256": sha((MATDIR / "replay_summary.json").read_bytes()),
        "original_pdf_sha256": OLD_HASH,
        "retained_circuit_stream_sha256": schematic_sha,
        "revised_pdf_sha256": sha(pdf),
        "figure_size_pt": [W, H],
        "placed_width_cm": 14,
        "compared_replay": str(args.replay) if args.replay else None,
        "main_view_s": MAIN[2:],
        "magnified_view_s": DETAIL[2:],
        "axes": {
            row: {
                key: {
                    "left": col[0],
                    "right": col[1],
                    "onset_x_pt": xy(col, row, 0, 0)[0],
                    "pt_per_s": (col[1] - col[0]) / (col[3] - col[2]),
                }
                for key, col in [("main", MAIN), ("detail", DETAIL)]
            }
            for row in ROWS
        },
        "numeric": numeric,
        "pdf_verified_trace_vertices": records,
    }
    (args.output / "FIGURE_6_4_NUMERICAL_VERIFICATION.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
