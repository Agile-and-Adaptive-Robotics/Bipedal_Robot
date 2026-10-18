import sys, statistics as st
folder = sys.argv[1]
rows = []
with open(folder + "/Rhythm Generator.txt", encoding="utf-8", errors="replace") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) == len(hdr):
            try: rows.append([float(v) for v in p])
            except ValueError: pass
while rows and all(v == 0.0 for v in rows[-1][2:]): rows.pop()
t = [r[1] for r in rows]
print(folder, "rows:", len(rows), "t_end:", t[-1] if t else None)
for col in ("L RG flx", "R RG flx", "L RG ext", "R RG ext"):
    ci = hdr.index(col)
    v = [r[ci]*1000 for r in rows]
    w = [vi for ti, vi in zip(t, v) if ti >= 0.5]  # post-transient window
    line = f"  {col:9s} post0.5s min {min(w):7.2f} max {max(w):7.2f} mean {st.mean(w):7.2f}"
    for thr in (-55.0, -56.0, -57.0, -58.0):
        up = sum(1 for a, b in zip(v, v[1:]) if a <= thr < b and True)
        # count only crossings at t>=0.3
        up2 = sum(1 for i in range(1, len(v)) if v[i-1] <= thr*1e-3 and False)
        n = sum(1 for i in range(1, len(v)) if t[i] >= 0.3 and v[i-1] <= thr and v[i] > thr)
        line += f" | x{thr:.0f}: {n}"
    print(line)
