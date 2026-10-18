import json, glob, csv, collections

D = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\tmp\gate_20261016"

def near(v, t, tol=1e-6):
    return v is not None and abs(v - t) < tol

for pat in ("round8_w*.jsonl", "round6_w*.jsonl"):
    rows = []
    for f in sorted(glob.glob(D + "\\" + pat)):
        for ln in open(f, encoding="utf-8", errors="replace"):
            ln = ln.strip()
            if ln.startswith("{"):
                rows.append(json.loads(ln))
    print(f"\n===== {pat}: {len(rows)} trials =====")
    sent = collections.Counter()
    for r in rows:
        v = r.get("value")
        m = r.get("metrics") or {}
        if m.get("nan"): sent["nan"] += 1; continue
        if near(v, -25.0): sent["frozen -25"] += 1; continue
        if near(v, -320.0): sent["no-cycles -320"] += 1; continue
        if near(v, -60.0): sent["bil-gate -60"] += 1; continue
        sent["genuine"] += 1
    print("objective buckets:", dict(sent))
    # bilateral = >=3+3 cycles per report
    bil = []
    for r in rows:
        m = r.get("metrics") or {}
        cr, cl = m.get("k_n_cycles_r"), m.get("k_n_cycles_l")
        v = r.get("value")
        if v is None or m.get("nan"): continue
        if near(v, -25.0) or near(v, -320.0) or near(v, -60.0): continue
        if (cr or 0) >= 3 and (cl or 0) >= 3:
            bil.append((v, r.get("number"), cr, cl, m.get("k_duty_r"), m.get("k_duty_l"), r.get("params", {}).get("drive_scale")))
    bil.sort()
    print(f"bilateral (>=3+3 cyc) non-sentinel: {len(bil)} / {len(rows)}")
    for v, n, cr, cl, dr, dl, s in bil[:6]:
        print(f"  trial={n} obj={v:.5f} cyc={cr}/{cl} duty={None if dr is None else round(dr,3)}/{None if dl is None else round(dl,3)} dscale={None if s is None else round(s,3)}")
    atmin = [b for b in bil if b[2] == 3 or b[3] == 3]
    print(f"bilateral trials sitting AT the 3-cycle gate minimum: {len(atmin)}")

print("\n===== wean_trials.csv (robot) =====")
with open(D + r"\wean_trials.csv", newline="", encoding="utf-8", errors="replace") as fh:
    rd = csv.DictReader(fh)
    rows = list(rd)
print("total rows:", len(rows))
tags = collections.Counter(r["tag"].rsplit("_", 1)[0] for r in rows)
print("tags:", dict(tags))
def fnum(r, k):
    try: return float(r[k])
    except Exception: return None
for pref in sorted(tags):
    sub = [r for r in rows if r["tag"].startswith(pref)]
    objs = [(fnum(r, "obj"), r) for r in sub]
    dead = sum(1 for o, r in objs if o is not None and o <= 5.001)
    top = sorted([x for x in objs if x[0] is not None], key=lambda x: -x[0])[:2]
    tops = "; ".join(f"t{r['tag'].rsplit('_',1)[-1]} obj={o:.2f} cyc={r['cycL']}/{r['cycR']} surv={r['survive']}" for o, r in top)
    print(f"{pref}: {len(sub)} rows, dead(obj<=5): {dead}, best: {tops}")
