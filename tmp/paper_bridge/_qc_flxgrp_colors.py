from PIL import Image
import sys, io, collections
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
im = Image.open(r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge\_qc_flxgrp_render.png").convert("RGB")
w, h = im.size
print("size:", w, h)
cnt = collections.Counter(im.getdata())
# report the most common NON-grayscale colors (saturation = max-min > 60)
sat = [(c, n) for c, n in cnt.items() if max(c) - min(c) > 60 and n > 20]
sat.sort(key=lambda x: -x[1])
print("top saturated colors (RGB, count):")
for c, n in sat[:14]:
    print("  ", c, n)
# marker fill test: find blue-ish pixels; check a small bounding region for enclosed white
blues = [(x, y) for y in range(h) for x in range(w)
         if (lambda p: p[2] > 120 and p[2] - p[1] > 60 and p[2] - p[0] > 60)(im.getpixel((x, y)))]
print("blue-ish pixel count:", len(blues))
if blues:
    xs = [p[0] for p in blues]; ys = [p[1] for p in blues]
    print("blue bbox:", min(xs), max(xs), min(ys), max(ys))
    # sample a dense blue cluster center and look at its interior for white
    cc = collections.Counter(blues)
    cx, cy = cc.most_common(1)[0][0]
    print("densest blue pixel at", (cx, cy), "color", im.getpixel((cx, cy)))
    white = sum(1 for dx in range(-3, 4) for dy in range(-3, 4)
                if (lambda p: p[0] > 220 and p[1] > 220 and p[2] > 220)(im.getpixel(((cx+dx) % w, (cy+dy) % h))))
    print("white pixels within 7x7 of densest blue:", white, "(high => hollow/open markers)")
