"""Measure x-tick digit heights: axes line = longest central dark horizontal
run within the BOTTOM 30% of the image; glyphs = dark pixels 4-45 px below it."""
import numpy as np
from PIL import Image

def measure(path, dpi, print_in=6.2):
    im = np.asarray(Image.open(path).convert("L"))
    h, w = im.shape
    dark = im < 128
    core = dark[int(h*0.70):, int(w*0.2):int(w*0.8)]
    counts = core.sum(axis=1)
    local = int(np.argmax(counts))
    axis_row = int(h*0.70) + local
    strip = dark[axis_row+4:min(axis_row+4+45, h), :]
    rows_with = strip.any(axis=1)
    runs, run = [], 0
    for r in rows_with:
        if r:
            run += 1
        else:
            if run:
                runs.append(run)
            run = 0
    if run:
        runs.append(run)
    native_in = w / dpi
    scale = print_in / native_in
    if runs:
        med = float(np.median(runs))
        pt_native = med / 0.7 / dpi * 72
        return axis_row, runs, med, pt_native, pt_native*scale, scale
    return axis_row, runs, None, None, None, scale

for name, path, dpi in [
    ("cam_1", r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_1.png", 150),
    ("cam_2", r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_2.png", 150),
    ("cam_3", r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_3.png", 150),
    ("cam_4", r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_4.png", 150),
    ("route8", r"D:\Github\Bipedal_Robot\tmp\route_check_run8_20261014.png", 110),
    ("route8b", r"D:\Github\Bipedal_Robot\tmp\route_check_run8b_20261014.png", 110),
]:
    axis_row, runs, med, pt_n, pt_p, scale = measure(path, dpi)
    if med:
        print(f"{name}: bottom axis row={axis_row}, glyph runs px={sorted(runs)}, median={med:.0f}px "
              f"-> ~{pt_n:.1f}pt native -> ~{pt_p:.1f}pt at 6.2in print (scale {scale:.2f}x)")
    else:
        print(f"{name}: bottom axis row={axis_row}, no glyphs below (runs={runs})")
