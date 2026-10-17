"""Refine crop of band-cropped figure PDFs: ink-bbox -> new CropBox -> MediaBox=CropBox.

Renders the full MediaBox canvas, finds ink, intersects with the current
CropBox, applies the refined CropBox via gs, then normalizes with
-dUseCropBox so the output MediaBox equals the CropBox (tight page).
"""
import re
import subprocess
import sys
import os
from PIL import Image

MIK = r"D:\MiKTeX\miktex\bin\x64"
DPI = 150
PAD = 5.0  # PDF points of padding around ink


def boxes(pdf):
    out = subprocess.run(
        [os.path.join(MIK, "pdfinfo.exe"), "-box", pdf],
        capture_output=True, text=True,
    ).stdout
    mb = re.search(r"MediaBox:\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", out)
    cb = re.search(r"CropBox:\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", out)
    return (
        tuple(float(v) for v in mb.groups()),
        tuple(float(v) for v in cb.groups()) if cb else None,
    )


def ink_bbox(png):
    im = Image.open(png).convert("L")
    w, h = im.size
    px = im.load()
    xmin, ymin, xmax, ymax = w, h, -1, -1
    for y in range(h):
        for x in range(w):
            if px[x, y] < 245:
                if x < xmin:
                    xmin = x
                if x > xmax:
                    xmax = x
                if y < ymin:
                    ymin = y
                if y > ymax:
                    ymax = y
    return xmin, ymin, xmax, ymax, w, h


def process(pdf):
    (mbx0, mby0, mbx1, mby1), cb = boxes(pdf)
    png = pdf + "_r"
    subprocess.run(
        [os.path.join(MIK, "pdftoppm.exe"), "-png", "-r", str(DPI), pdf, png],
        check=True, capture_output=True,
    )
    png += "-1.png"
    x0, y0, x1, y1, w, h = ink_bbox(png)
    os.remove(png)
    if x1 < 0:
        print(pdf, ": EMPTY")
        return
    # px -> user space (MediaBox origin bottom-left; render origin top-left)
    sx = (mbx1 - mbx0) / w
    sy = (mby1 - mby0) / h
    left = mbx0 + x0 * sx - PAD
    right = mbx0 + x1 * sx + PAD
    top = mby1 - y0 * sy + PAD
    bottom = mby1 - y1 * sy - PAD
    cur = cb if cb else (mbx0, mby0, mbx1, mby1)
    left, bottom = max(left, cur[0]), max(bottom, cur[1])
    right, top = min(right, cur[2]), min(top, cur[3])
    box = f"[{left:.2f} {bottom:.2f} {right:.2f} {top:.2f}]"
    tmp = pdf + "_c.pdf"
    subprocess.run(
        [os.path.join(MIK, "mgs.exe"), "-q", "-o", tmp, "-sDEVICE=pdfwrite",
         "-c", f"[/CropBox {box} /PAGES pdfmark", "-f", pdf],
        check=True, capture_output=True,
    )
    out = pdf.replace(".pdf", "_fit.pdf")
    subprocess.run(
        [os.path.join(MIK, "mgs.exe"), "-q", "-o", out, "-sDEVICE=pdfwrite",
         "-dUseCropBox", tmp],
        check=True, capture_output=True,
    )
    os.remove(tmp)
    _, ncb = boxes(out)
    print(os.path.basename(out), "->", box, "| final box:", ncb)


for p in sys.argv[1:]:
    process(p)
