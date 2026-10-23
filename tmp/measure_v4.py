import os

from PIL import Image

os.chdir(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\renders")
im = Image.open("v4_p164.png").convert("L")
w, h = im.size
px = im.load()
minx, maxx = w, -1
miny, maxy = h, -1
for y in range(0, h, 2):
    for x in range(0, w, 2):
        if px[x, y] < 245:
            minx = min(minx, x)
            maxx = max(maxx, x)
            miny = min(miny, y)
            maxy = max(maxy, y)
print("v4_p164 ink x[", minx, ",", maxx, "] y[", miny, ",", maxy, "] of", w, "right-margin", w - maxx)
