"""Measure tick-digit pixel heights: find the axes bottom line (longest dark
horizontal run), then measure connected-component heights of dark pixels in the
band 4-60 px below it (x tick labels). Also left-side band for y ticks."""
import numpy as np
from PIL import Image

def digit_heights(path, band=60):
    im = np.asarray(Image.open(path).convert("L"))
    h, w = im.shape
    dark = im < 128
    # longest horizontal dark run per row -> the axis line
    row_counts = dark.sum(axis=1)
    # restrict to central columns to avoid border artifacts
    core = dark[:, int(w*0.2):int(w*0.8)]
    core_counts = core.sum(axis=1)
    axis_row = int(np.argmax(core_counts))
    # search downward for glyph rows
    strip = dark[axis_row+4:min(axis_row+4+band, h), :]
    rows_with = strip.any(axis=1)
    heights, run = [], 0
    for r in rows_with:
        if r:
            run += 1
        else:
            if run:
                heights.append(run)
            run = 0
    if run:
        heights.append(run)
    return axis_row, heights

for name, path, dpi, native_w_in in [
    ("cam_1", r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_1.png", 150, None),
    ("cam_2", r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_2.png", 150, None),
    ("cam_3", r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_3.png", 150, None),
    ("cam_4", r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_4.png", 150, None),
    ("route8", r"D:\Github\Bipedal_Robot\tmp\route_check_run8_20261014.png", 110, None),
    ("route8b", r"D:\Github\Bipedal_Robot\tmp\route_check_run8b_20261014.png", 110, None),
]:
    im = Image.open(path)
    w_px = im.size[0]
    native_in = w_px / dpi
    scale = 6.2 / native_in  # print at 6.2 in
    axis_row, hs = digit_heights(path)
    if hs:
        med = float(np.median(hs))
        # digit cap height ~0.7 of point size: pt_native = med/0.7/dpi*72
        pt_native = med / 0.7 / dpi * 72
        print(f"{name}: axis_row={axis_row} glyph-run heights px={sorted(hs)} median={med:.0f}px "
              f"-> ~{pt_native:.1f}pt native, ~{pt_native*scale:.1f}pt at 6.2in print (scale {scale:.2f}x)")
    else:
        print(f"{name}: axis_row={axis_row} NO glyphs found below axis line")
