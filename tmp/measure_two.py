import os

from PIL import Image

os.chdir(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\renders")
for f in ("base_p155.png", "base_p204.png"):
    im = Image.open(f).convert("L")
    w, h = im.size
    px = im.load()
    minx, miny, maxx, maxy = w, h, -1, -1
    for y in range(0, h, 2):
        for x in range(0, w, 2):
            if px[x, y] < 245:
                minx = min(minx, x)
                maxx = max(maxx, x)
                miny = min(miny, y)
                maxy = max(maxy, y)
    print(f, "page", w, "x", h, "ink x[", minx, ",", maxx, "] y[", miny, ",", maxy, "]")
