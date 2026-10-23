import glob
import io
import sys

tree = sys.argv[1]
ABS = "D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Figures/"
for p in glob.glob(tree + "/**/*.tex", recursive=True):
    with io.open(p, "r", encoding="utf-8-sig", newline="") as f:
        t = f.read()
    if "../Figures/" in t:
        t = t.replace("../Figures/", ABS)
        with io.open(p, "w", encoding="utf-8", newline="") as f:
            f.write(t)
        print("fixed", p)
