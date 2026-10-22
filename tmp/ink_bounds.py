"""Ink-bounds check: for each rendered page PNG, find the non-white bounding box
and flag pages whose ink extends unusually close to the page edge (margin spill)
compared with control pages from untouched chapters."""
import glob
import os

from PIL import Image

REND = os.path.join(
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation"
    r"\Notes\advisor_revision_20261022\renders"
)

# control pages rendered earlier from the proposal build (untouched chapters):
# p101 (ch4 body), p220 (appendix) were rendered in the first pass?  We only
# have the first-pass set (wrong offsets) + v2 set.  Use first-pass pages from
# untouched regions as controls: those were PDF pages 14,15,17,32,33,35,37,38,
# 42,43,129,132,133,136,137,142,143,145,147,150,151,168,172,185 (old offsets
# => they landed in untouched material for many).
controls = [os.path.join(REND, f"p{n}.png") for n in (14, 15, 32, 33, 35, 37, 38, 42, 43, 129)]
targets = sorted(glob.glob(os.path.join(REND, "v2_p*.png")))


def ink_bbox(path, thresh=245):
    im = Image.open(path).convert("L")
    w, h = im.size
    px = im.load()
    minx, miny, maxx, maxy = w, h, -1, -1
    for y in range(0, h, 2):
        for x in range(0, w, 2):
            if px[x, y] < thresh:
                if x < minx:
                    minx = x
                if x > maxx:
                    maxx = x
                if y < miny:
                    miny = y
                if y > maxy:
                    maxy = y
    return w, h, minx, miny, maxx, maxy


print("CONTROLS (untouched pages):")
ctrl_right = []
for p in controls:
    if not os.path.exists(p):
        continue
    w, h, x0, y0, x1, y1 = ink_bbox(p)
    ctrl_right.append((x1, w))
    print(f"  {os.path.basename(p):12s} page={w}x{h} ink x[{x0},{x1}] y[{y0},{y1}] right-margin={w-x1}px")

print("\nTARGETS (revised pages):")
for p in targets:
    w, h, x0, y0, x1, y1 = ink_bbox(p)
    rm = w - x1
    flag = ""
    if ctrl_right:
        import statistics
        mu = statistics.mean(c[0] for c in ctrl_right)
        sd = statistics.pstdev(c[0] for c in ctrl_right)
        if x1 > mu + 3 * sd:
            flag = "  <<< RIGHT EDGE BEYOND CONTROLS"
    if rm < 20:
        flag += "  <<< VERY CLOSE TO PAGE EDGE"
    print(f"  {os.path.basename(p):14s} page={w}x{h} ink x[{x0},{x1}] y[{y0},{y1}] right-margin={rm}px{flag}")
