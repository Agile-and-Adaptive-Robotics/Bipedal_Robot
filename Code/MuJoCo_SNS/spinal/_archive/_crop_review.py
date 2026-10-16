"""Crop a rendered figure into native-resolution review regions.

Usage: python _crop_review.py <png> <outprefix> x0,y0,x1,y1 [...]
Coordinates are FRACTIONS of the image (x right, y DOWN), so crops work
independently of the tight-bbox output size. Writes <outprefix>_<n>.png.
"""
import sys

from PIL import Image


def main():
    src = Image.open(sys.argv[1])
    W, H = src.size
    print(f"{sys.argv[1]}: {W}x{H}")
    out = sys.argv[2]
    i = 3
    n = 0
    while i < len(sys.argv):
        x0, y0, x1, y1 = (float(v) for v in sys.argv[i].split(","))
        box = (int(x0 * W), int(y0 * H), int(x1 * W), int(y1 * H))
        p = f"{out}_{n}.png"
        src.crop(box).save(p)
        print(f"{p}  box={box}")
        n += 1
        i += 1


if __name__ == "__main__":
    main()
