"""Print per-line text layout with y-coordinates from pdftotext -bbox output."""
import re
import sys

def analyze(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    words = re.findall(
        r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>',
        txt,
    )
    if not words:
        print(f"{path}: NO WORDS")
        return
    # group words into lines by yMin (rows within 3pt)
    lines = {}
    for x0, y0, x1, y1, w in words:
        y0f = float(y0)
        key = None
        for k in lines:
            if abs(k - y0f) < 3.5:
                key = k
                break
        if key is None:
            key = y0f
            lines[key] = []
        lines[key].append((float(x0), float(x1), w))
    print(f"===== {path.split(chr(92))[-1]}")
    for y in sorted(lines):
        ws = sorted(lines[y])
        text = " ".join(w for _, _, w in ws)
        xmin = min(x0 for x0, _, _ in ws)
        xmax = max(x1 for _, x1, _ in ws)
        if len(text) > 90:
            text = text[:90] + "..."
        print(f"y={y:7.1f}-{y:7.1f} x=[{xmin:6.1f},{xmax:6.1f}] {text}")

for p in sys.argv[1:]:
    analyze(p)
