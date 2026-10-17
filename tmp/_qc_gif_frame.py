"""QC helper: dump GIF frame count/size + first frame to PNG; report PNG sizes/DPI."""
import os
from PIL import Image

root = r"D:\Github\Bipedal_Robot"
targets = [
    r"Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_1.png",
    r"Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_2.png",
    r"Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_3.png",
    r"Code\Matlab\Mesh_Optimization\Results\Vas_Pam_20mm_cam_winner_4.png",
    r"tmp\route_check_run8_20261014.png",
    r"tmp\route_check_run8b_20261014.png",
]

for rel in targets:
    p = os.path.join(root, rel)
    im = Image.open(p)
    print(f"{rel}: {im.size[0]}x{im.size[1]} px, mode={im.mode}, dpi={im.info.get('dpi', 'none')}")

gif = os.path.join(root, r"Code\Matlab\Mesh_Optimization\Results\Knee_Extensor_20mm.gif")
g = Image.open(gif)
n = getattr(g, "n_frames", 1)
print(f"Knee_Extensor_20mm.gif: {g.size[0]}x{g.size[1]} px, frames={n}, mode={g.mode}")
g.seek(0)
first = g.convert("RGB")
out = os.path.join(root, r"tmp\_qc_gif_first_frame.png")
first.save(out)
print(f"first frame saved -> {out}")
