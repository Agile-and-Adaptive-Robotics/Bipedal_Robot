"""Crop a figure region out of a PDF page as vector PDF (pdfwrite + CropBox),
with a PNG preview for verification. Usage:
  python _crop.py <pdf> <page> <x0f> <y0f> <x1f> <y1f> <out_pdf>
Fractions are of page width/height, y measured from TOP.
Step 1 (-p flag): produce a tight ink-trim analysis of the cropped preview.
"""
import re
import subprocess
import sys
import tempfile

from PIL import Image

MGS = r"D:\MiKTeX\miktex\bin\x64\mgs.exe"


def pdf_page_size(pdf, page):
    """True MediaBox size in points, via a calibration render at 100 DPI
    (rendering uses the full un-cropped page)."""
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        out = f.name
    subprocess.run([MGS, "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m",
                    "-r100", f"-dFirstPage={page}", f"-dLastPage={page}",
                    f"-sOutputFile={out}", pdf],
                   check=True, capture_output=True)
    im = Image.open(out)
    w_pt = im.size[0] / 100 * 72
    h_pt = im.size[1] / 100 * 72
    return [0.0, 0.0, w_pt, h_pt]


def crop_pdf(pdf, page, fracs, out_pdf):
    x0, y0, x1, y1 = pdf_page_size(pdf, page)
    W, H = x1 - x0, y1 - y0
    # y from top -> PDF y from bottom
    bx0 = x0 + W * fracs[0]
    bx1 = x0 + W * fracs[2]
    by1 = y1 - H * fracs[1]   # top edge (fraction from top)
    by0 = y1 - H * fracs[3]   # bottom edge
    cmd = [MGS, "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
           f"-dFirstPage={page}", f"-dLastPage={page}",
           f"-o{out_pdf}",
           "-c", f"[/CropBox [{bx0:.2f} {by0:.2f} {bx1:.2f} {by1:.2f}]"
                 " /PAGES pdfmark",
           "-f", pdf]
    subprocess.run(cmd, check=True, capture_output=True)
    return out_pdf


def render(pdf, dpi, out_png):
    subprocess.run([MGS, "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m",
                    "-dUseCropBox", f"-r{dpi}", f"-sOutputFile={out_png}",
                    pdf],
                   check=True, capture_output=True)


def ink_bbox(png, thresh=245, pad=6):
    im = Image.open(png).convert("L")
    px = im.load()
    w, h = im.size
    minx, miny, maxx, maxy = w, h, -1, -1
    for y in range(h):
        for x in range(w):
            if px[x, y] < thresh:
                if x < minx:
                    minx = x
                if x > maxx:
                    maxx = x
                if y < miny:
                    miny = y
                if y > maxy:
                    maxy = y
    if maxx < 0:
        return None
    return (max(0, minx - pad), max(0, miny - pad),
            min(w, maxx + pad), min(h, maxy + pad))


if __name__ == "__main__":
    args = sys.argv[1:]
    mode = args[0]
    if mode == "crop":
        pdf, page = args[1], int(args[2])
        fracs = [float(v) for v in args[3:7]]
        out_pdf, out_png = args[7], args[8]
        crop_pdf(pdf, page, fracs, out_pdf)
        render(out_pdf, 150, out_png)
        im = Image.open(out_png)
        print("cropped", out_pdf, "preview", im.size)
    elif mode == "trim":
        # report ink bbox of a preview png as fractions of its size
        png = args[1]
        bb = ink_bbox(png)
        im = Image.open(png)
        w, h = im.size
        print(f"ink bbox px {bb} of {w}x{h} -> "
              f"fractions x0={bb[0]/w:.4f} y0={bb[1]/h:.4f} "
              f"x1={bb[2]/w:.4f} y1={bb[3]/h:.4f}")
    elif mode == "png":
        # high-res PNG from an already-cropped PDF
        pdf, dpi, out = args[1], int(args[2]), args[3]
        render(pdf, dpi, out)
        print("wrote", out, Image.open(out).size)
