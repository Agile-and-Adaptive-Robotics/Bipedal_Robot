"""Supervisor M2 audit — independent figure QA (read-only)."""
from PIL import Image
import numpy as np
import re
import os

os.chdir(os.path.join("D:", os.sep, "GitHub", "Bipedal_Robot",
                      "Documentation", "Reports and Papers", "Dissertation",
                      "Figures"))

figs = [
    "30-results/spiking_mirror_rhythm",
    "30-results/spiking_mirror_ground_gait",
    "30-results/simulink_knee_spiking_twin",
    "93-AppendixB/animatlab_spiking_conversion",
]
tol = {"indigo#0000FF": (0, 0, 255), "orange#FFB14E": (255, 177, 78),
       "pink#EA5F94": (234, 95, 148)}

for f in figs:
    png, pdf, alt = f + ".png", f + ".pdf", f + "_alt.txt"
    im = Image.open(png)
    a = np.asarray(im.convert("RGB"))
    h, w, _ = a.shape
    nonwhite = (a.astype(int).sum(axis=2) < 720)
    rows = np.where(nonwhite.any(axis=1))[0]
    cols = np.where(nonwhite.any(axis=0))[0]
    margins = (int(rows[0]), int(h - 1 - rows[-1]), int(cols[0]),
               int(w - 1 - cols[-1]))  # t/b/l/r
    flat = a.reshape(-1, 3).astype(int)
    present = {}
    for name, c in tol.items():
        d = np.abs(flat - np.array(c)).max(axis=1)
        present[name] = int((d <= 2).sum())
    raw = open(pdf, "rb").read()
    fonts = sorted(set(x.decode() for x in
                       re.findall(rb"/BaseFont\s*/([A-Za-z0-9+\-]+)", raw)))
    italic = sorted(set(float(x) for x in
                        re.findall(rb"/ItalicAngle\s+(-?[\d.]+)", raw)))
    asz = os.path.getsize(alt) if os.path.exists(alt) else -1
    print(os.path.basename(png))
    print("  size px %dx%d = %.2fx%.2f in @300dpi | margins t/b/l/r %s | alt.txt %d bytes"
          % (w, h, w / 300.0, h / 300.0, margins, asz))
    print("  tol-color pixels %s" % present)
    print("  fonts %s | ItalicAngle %s" % (fonts, italic))
