from PIL import Image
import sys, io, collections
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
for name in ["_qc_FlxPin_group.png", "_qc_ExtPin_group.png", "_qc_Ext10mm_52cm.png"]:
    im = Image.open(r"D:\GitHub\Bipedal_Robot\tmp\paper_bridge" + "\\" + name).convert("RGB")
    w, h = im.size
    cnt = collections.Counter(im.getdata())
    sat = [(c, n) for c, n in cnt.items() if max(c) - min(c) > 60 and n > 15]
    sat.sort(key=lambda x: -x[1])
    print(f"\n{name} ({w}x{h}) top saturated colors:")
    for c, n in sat[:10]:
        print("   RGB", c, "count", n)
    # indigo test: dark blue with significant red component (75,0,130)-ish
    indigo = sum(n for c, n in cnt.items() if c[2] > 90 and c[2] - c[1] > 40 and c[2] - c[0] > 30 and c[0] > 25 and c[0] < 140)
    pureblue = sum(n for c, n in cnt.items() if c[2] > 180 and c[0] < 60 and c[1] < 60)
    gold = sum(n for c, n in cnt.items() if c[0] > 200 and 140 < c[1] < 215 and c[2] < 110)
    gray = sum(n for c, n in cnt.items() if abs(c[0]-102) < 12 and abs(c[1]-102) < 12 and abs(c[2]-102) < 12)
    print(f"   indigo-ish px: {indigo} | pure-blue px: {pureblue} | gold-ish px: {gold} | gray(0.4) px: {gray}")
