"""Render a PDF page range at low DPI into a labeled contact sheet grid."""
import subprocess
import sys
import tempfile
from PIL import Image, ImageDraw

MGS = r"D:\MiKTeX\miktex\bin\x64\mgs.exe"


def render(pdf, page, dpi, out):
    subprocess.run([MGS, "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=png16m",
                    f"-r{dpi}", f"-dFirstPage={page}", f"-dLastPage={page}",
                    f"-sOutputFile={out}", pdf],
                   check=True, capture_output=True)


pdf = sys.argv[1]
p0, p1 = int(sys.argv[2]), int(sys.argv[3])
cols = int(sys.argv[4]) if len(sys.argv) > 4 else 4
sheet = sys.argv[5] if len(sys.argv) > 5 else "sheet.png"
dpi = int(sys.argv[6]) if len(sys.argv) > 6 else 60

thumbs = []
for p in range(p0, p1 + 1):
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        out = f.name
    render(pdf, p, dpi, out)
    im = Image.open(out).convert("RGB")
    d = ImageDraw.Draw(im)
    d.rectangle([0, 0, 78, 26], fill="yellow")
    d.text((6, 6), f"p{p}", fill="black")
    thumbs.append(im)

w, h = thumbs[0].size
rows = (len(thumbs) + cols - 1) // cols
grid = Image.new("RGB", (w * cols, h * rows), "white")
for i, im in enumerate(thumbs):
    grid.paste(im, ((i % cols) * w, (i // cols) * h))
grid.save(sheet)
print("wrote", sheet, grid.size)
