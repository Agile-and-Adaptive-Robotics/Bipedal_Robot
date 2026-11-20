"""Draw the signed flexor data on the native extensor publication panel."""

import argparse
import json
from math import sqrt
from pathlib import Path

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas


PAGE = (422, 338)
LEFT, RIGHT, BOTTOM, TOP = 46.5, 418.5, 39.125, 332.375
XLIM = (-125, 15)
INK = (0.1294, 0.1294, 0.1294)
BLACK = (0, 0, 0)
MAGENTA = (0.804, 0.204, 0.710)
GOLD = (1, 0.843, 0)
BLUE = (0, 0, 1)


def scale(value, limits, start, end):
    return start + (value - limits[0]) * (end - start) / (limits[1] - limits[0])


def draw_axes(pdf, x, y):
    pdf.setStrokeColorRGB(*INK)
    pdf.setFillColorRGB(*INK)
    pdf.setLineWidth(2)
    pdf.setFont("ArialBold", 12)
    pdf.line(LEFT, BOTTOM, RIGHT, BOTTOM)
    pdf.line(LEFT, BOTTOM, LEFT, TOP)
    for value in range(XLIM[0], XLIM[1] + 1, 5):
        major = value % 30 == 0
        pdf.line(x(value), BOTTOM, x(value), BOTTOM + (9.3 if major else 4.65))
        if major:
            pdf.drawCentredString(x(value), BOTTOM - 15.65, str(value))
    for value in range(-104, 5, 2):
        major = value % 10 == 0
        pdf.line(LEFT, y(value), LEFT + (9.3 if major else 4.65), y(value))
        if major:
            pdf.drawRightString(LEFT - 4.5, y(value) - 4.12, str(value))
    pdf.saveState()
    pdf.translate(9.5, (TOP + BOTTOM) / 2)
    pdf.rotate(90)
    pdf.setFont("ArialBold", 10)
    pdf.drawCentredString(0, 0, "Torque, N m")
    pdf.restoreState()
    # The original letters are theta, subscript k, comma, and degree sign.
    cx = (LEFT + RIGHT) / 2
    pdf.setFont("ArialBoldItalic", 12.75)
    pdf.drawString(cx - 13.8, 8.0, "θ")
    pdf.setFont("ArialBold", 10.5)
    pdf.drawString(cx - 5.9, 2.5, "k")
    pdf.setFont("ArialBold", 12.75)
    pdf.drawString(cx + 1, 8.0, ", °")


def stroke_curve(pdf, points, role):
    pdf.saveState()
    if role == "human":
        pdf.setStrokeColorRGB(*BLACK)
        pdf.setLineWidth(3)
        pdf.setDash([3, 3], 0)
    else:
        pdf.setStrokeColorRGB(*MAGENTA)
        pdf.setLineWidth(2.5)
    pen = pdf.beginPath()
    pen.moveTo(*points[0])
    for pt in points[1:]:
        pen.lineTo(*pt)
    pdf.drawPath(pen, stroke=1, fill=0)
    pdf.restoreState()


def marker(pdf, x, y, area, color, alpha):
    pdf.saveState()
    pdf.setFillColorRGB(*color)
    pdf.setFillAlpha(alpha)
    pdf.circle(x, y, sqrt(area) / 2, stroke=0, fill=1)
    pdf.restoreState()


def legend(pdf, data):
    # This corner is below the plotted data for the full width of the key.
    x0, y0, width, height = 64.50, 55.375, 194.25, 48.75
    pdf.setLineWidth(2)
    pdf.setFillColorRGB(1, 1, 1)
    pdf.setStrokeColorRGB(*INK)
    pdf.rect(x0, y0, width, height, fill=1, stroke=1)
    pdf.setFont("ArialBold", 9.75)
    pdf.setFillColorRGB(*BLACK)
    entries = (
        ("design", "Two-BPA design, 620 kPa"),
        ("human", "Biceps femoris short-head target"),
        ("prediction", "Prediction at measured pressure"),
        ("measured", "Measured"),
    )
    for index, (role, label) in enumerate(entries):
        baseline = y0 + height - 10.45 - index * 11.535
        cy = baseline + 3.38
        if role in ("design", "human"):
            stroke_curve(pdf, [(x0 + 3, cy), (x0 + 33, cy)], role)
        else:
            rgb = BLUE if role == "measured" else GOLD
            marker(pdf, x0 + 18, cy, 36, rgb, data[role]["alpha"])
        assert pdfmetrics.stringWidth(label, "ArialBold", 9.75) < width - 36
        pdf.drawString(x0 + 35.25, baseline, label)


def draw(data_path, output_path, font_dir):
    source = json.loads(data_path.read_text(encoding="utf-8"))
    data = {series["name"]: series for series in source["series"]}
    assert sorted(len(s["x"]) for s in data.values()) == [15, 15, 100, 100]
    assert source["ylim"] == [-105, 5]
    pdfmetrics.registerFont(TTFont("ArialBold", font_dir / "arialbd.ttf"))
    pdfmetrics.registerFont(TTFont("ArialBoldItalic", font_dir / "arialbi.ttf"))
    pdf = canvas.Canvas(str(output_path), pagesize=PAGE, invariant=1,
                        initialFontName="ArialBold", initialFontSize=12)
    pdf.setTitle("Two-BPA flexor validation")
    pdf.setCreator("Bipedal_Robot physical validation figure generator")
    x = lambda angle: scale(angle, XLIM, LEFT, RIGHT)
    y = lambda torque: scale(torque, source["ylim"], BOTTOM, TOP)
    draw_axes(pdf, x, y)
    for role in ("design", "human"):
        s = data[role]
        stroke_curve(pdf, [(x(a), y(t)) for a, t in zip(s["x"], s["y"])], role)
    for role, rgb in (("prediction", GOLD), ("measured", BLUE)):
        s = data[role]
        for a, t in zip(s["x"], s["y"]):
            marker(pdf, x(a), y(t), s["area"], rgb, s["alpha"])
    legend(pdf, data)
    pdf.showPage()
    pdf.save()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("--data", type=Path, default=Path(__file__).with_name("flexor_data.json"))
    parser.add_argument("--font-dir", type=Path, default=Path("C:/Windows/Fonts"))
    args = parser.parse_args()
    draw(args.data, args.output_path, args.font_dir)
