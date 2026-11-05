"""make_gifs.py — assemble looping GIFs from rendered frame sets.

The defense computer embeds the GIF files directly (no external links).
Each GIF loops forever, is quantized to a 256-color adaptive palette, and
is downscaled to 640 px to keep the deck file size sane.
"""
from pathlib import Path
from PIL import Image

HERE = Path(__file__).resolve().parents[1]
FPS = 15
GIF_W = 640

JOBS = [
    ("visual1", "cam_extensor_loop.gif"),
    ("visual2", "reverse_pulley_loop.gif"),
    ("visual3", "knee_frames_loop.gif"),
]

for folder, out_name in JOBS:
    fdir = HERE / "build" / folder / "frames"
    frames = sorted(fdir.glob("f*.png"))
    if not frames:
        print(f"[{folder}] NO FRAMES FOUND — skipped")
        continue
    ims = []
    for f in frames:
        im = Image.open(f).convert("RGB")
        if im.width != GIF_W:
            im = im.resize((GIF_W, GIF_W), Image.LANCZOS)
        ims.append(im)
    # quantize first frame fully, others fast (palette reuse)
    pal = ims[0].quantize(colors=256, method=Image.MEDIANCUT)
    qims = [im.quantize(palette=pal, dither=Image.FLOYDSTEINBERG) for im in ims]
    out = HERE / "build" / folder / out_name
    qims[0].save(out, save_all=True, append_images=qims[1:], duration=int(1000 / FPS),
                 loop=0, optimize=False)
    print(f"[{folder}] {len(qims)} frames -> {out.name}  {out.stat().st_size/1e6:.1f} MB")

print("GIFS_DONE")
