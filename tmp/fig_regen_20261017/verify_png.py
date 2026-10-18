"""Verify the saved render: dimensions, dpi metadata, white background."""
from pathlib import Path
import numpy as np
from PIL import Image

p = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\30-torque\opensim_robot_isometric.png")
im = Image.open(p)
print("size:", im.size, "dpi:", im.info.get("dpi"), "mode:", im.mode)
a = np.asarray(im.convert("RGB")).astype(int)
h, w, _ = a.shape
corners = dict(
    tl=a[5, 5].tolist(), tr=a[5, w - 6].tolist(),
    bl=a[h - 6, 5].tolist(), br=a[h - 6, w - 6].tolist(),
    top_mid=a[5, w // 2].tolist(),
)
print("corner pixels (want ~[255,255,255]):", corners)
near_white = (a > 246).all(axis=2)
print(f"near-white fraction: {near_white.mean():.3f}")
red = (a[:, :, 0] > 120) & (a[:, :, 1] < 90) & (a[:, :, 2] < 90)
print(f"red muscle-path pixel fraction: {red.mean():.4f}  (count={int(red.sum())})")
# any pure-black pixels (garbled-render symptom)?
black = (a < 10).all(axis=2)
print(f"near-black fraction: {black.mean():.5f}")
