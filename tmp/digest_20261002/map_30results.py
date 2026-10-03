t = open(r'C:\Users\Ben Bolen\AppData\Local\Temp\o808_30results.tex', encoding='utf8').read()
lines = t.splitlines()
anchors = ["sec:ongoing", "re-identified", "The adopted values", "far advanced",
           "produced usable", "ZCODE 2026-09-26", "ZCODE 2026-10-02", "\\section", "VU 0.67"]
for i, ln in enumerate(lines):
    for a in anchors:
        if a in ln:
            print(f"{i+1:5d}  [{a}]  {ln[:110]}")
            break
print("total lines:", len(lines))
