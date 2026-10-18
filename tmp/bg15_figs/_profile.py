"""Row/column ink profile of a rendered page to find content bands precisely.
Prints runs of rows whose ink fraction is below a threshold (gaps), so figure
blocks can be delimited from header/caption/body text.
"""
import sys

import numpy as np
from PIL import Image

png = sys.argv[1]
axis = sys.argv[2] if len(sys.argv) > 2 else "rows"
thr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.002
r0 = int(sys.argv[4]) if len(sys.argv) > 4 else None
r1 = int(sys.argv[5]) if len(sys.argv) > 5 else None

im = np.asarray(Image.open(png).convert("L"), dtype=np.uint8)
ink = (im < 245).astype(np.float32)
if r0 is not None:
    ink[:r0, :] = 0
    if r1 is not None:
        ink[r1 + 1:, :] = 0
if axis == "rows":
    prof = ink.mean(axis=1)
    n = len(prof)
else:
    prof = ink.mean(axis=0)
    n = len(prof)

# print contiguous bands of ink (>thr) as start-end fractions
bands = []
start = None
for i, v in enumerate(prof):
    if v > thr and start is None:
        start = i
    elif v <= thr and start is not None:
        bands.append((start, i - 1))
        start = None
if start is not None:
    bands.append((start, n - 1))

# merge bands separated by <3 units
merged = []
for b in bands:
    if merged and b[0] - merged[-1][1] <= 3:
        merged[-1] = (merged[-1][0], b[1])
    else:
        merged.append(list(b))

print(f"{png}: {n} px, {len(merged)} bands")
for a, b in merged:
    print(f"  {a}-{b} px  ({a/n:.4f}-{b/n:.4f})")
