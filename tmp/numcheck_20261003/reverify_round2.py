import scipy.io as sio
import numpy as np
import csv, json, os

R = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"

def load(f):
    return sio.loadmat(R + "\\" + f, squeeze_me=True, struct_as_record=False)

print("== 1. Extensor design of record (Vas_Pam_20mm_Result.mat) ==")
d = load("Vas_Pam_20mm_Result.mat")
tm = np.asarray(d["torqueMargin"], dtype=float).ravel()
ang = np.asarray(d["ctx"].humanAngleD, dtype=float).ravel()
tz = np.asarray(d["predBest"].TorqueZ, dtype=float).ravel()
hu = np.asarray(d["ctx"].humanTorque, dtype=float).ravel()
print("fBest=%.5f  minMargin=%.6f at 1-based idx %d -> angle %.4f deg" %
      (float(d["fBest"]), tm.min(), tm.argmin()+1, ang[tm.argmax()] if False else ang[tm.argmin()]))
print("TorqueZ span %.4f..%.4f  humanTorque span %.4f..%.4f  shortfalls below human: %d" %
      (tz.min(), tz.max(), hu.min(), hu.max(), int((tm < 0).sum())))
print("ctx: Dia=%s BPAcount=%s P=%s KMAX=%s" %
      (d["ctx"].Dia, d["ctx"].BPAcount, d["ctx"].targetPressure, d["ctx"].KMAX))

def find_field(obj, name, path="root", depth=0):
    """Recursively find a field by name in nested mat structs; return (path, value) list."""
    hits = []
    if depth > 4:
        return hits
    if isinstance(obj, sio.matlab.mio5_params.mat_struct):
        for f in obj._fieldnames:
            v = getattr(obj, f)
            if f == name:
                hits.append((path + "." + f, v))
            hits += find_field(v, name, path + "." + f, depth + 1)
    elif isinstance(obj, np.ndarray) and obj.dtype == object:
        for i, v in enumerate(obj.ravel()[:3]):
            hits += find_field(v, name, "%s[%d]" % (path, i), depth + 1)
    return hits

for f in ["XiUsed", "removedText"]:
    hits = find_field(d, f)
    for pth, v in hits:
        print("Vas %s at %s = %s" % (f, pth, v if not hasattr(v, "ravel") else np.asarray(v).ravel()))
print("keys with 'removed':", [k for k in d if "remove" in k.lower()])
if "removedText" in d:
    print("removedText =", d["removedText"])

print("\n== 2. 09-20 snapshot ==")
d2 = load("Vas_Pam_20mm_Result_20260920_1519.mat")
print("fBest=%.5f" % float(d2["fBest"]))

print("\n== 3. 09-25 dated mat (Appendix C) ==")
d3 = load("Vas_Pam_20mm_Result_20260925.mat")
xb = np.asarray(d3["xBest"], dtype=float).ravel()
ea = np.asarray(d3["routeInfo"].eliminatedAngleD, dtype=float).ravel()
eo = np.asarray(d3["routeInfo"].eliminationOrder, dtype=float).ravel().astype(int)
print("xBest =", xb)
print("p1_mm =", xb[0:3]*1000, " pEnd_mm =", xb[3:6]*1000, " rest_mm = %.1f tendon_mm = %.1f" % (xb[6]*1000, xb[7]*1000))
print("grid: n=%d %.2f..%.2f spacing %.4f" % (ang.size, ang.min(), ang.max(), ang[1]-ang[0]))
print("eliminatedAngleD =", ea, " order =", eo)
order = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
for i, e in enumerate(ea):
    if not np.isnan(e):
        panel = e - (ang[1]-ang[0])
        print("  %s released at %+.2f (panel before: %+.2f)" % (order[i], e, panel))
print("keys with 'removed':", [k for k in d3 if "remove" in k.lower()])
for f in ["XiUsed", "removedText", "removed"]:
    for pth, v in find_field(d3, f):
        print("Vas0925 %s at %s = %s" % (f, pth, v if not hasattr(v, "ravel") else np.asarray(v).ravel()))
pl0 = np.asarray(d3["predBest"].pathLength0, dtype=float).ravel()
i0 = int(np.argmin(np.abs(ang)))
print("pathLength0 at theta=%.2f: %.4f m" % (ang[i0], pl0[i0]))

print("\n== 4. Flexor (Bifemsh_20mm_Result.mat) ==")
d4 = load("Bifemsh_20mm_Result.mat")
tmf = np.asarray(d4["torqueMarginFraction"], dtype=float).ravel()
angf = np.asarray(d4["ctx"].humanAngleD, dtype=float).ravel()
tzf = np.asarray(d4["predBest"].TorqueZ, dtype=float).ravel()
huf = np.asarray(d4["ctx"].humanTorque, dtype=float).ravel()
print("fBest=%.4f Jseed=%.2f minMarginFraction=%.6f at %.2f deg" %
      (float(d4["fBest"]), float(d4["Jseed"]), tmf.min(), angf[tmf.argmin()]))
print("TorqueZ span %.4f..%.4f  humanTorque span %.4f..%.4f" % (tzf.min(), tzf.max(), huf.min(), huf.max()))
for pth, v in find_field(d4, "XiUsed"):
    print("Bifemsh XiUsed at %s = %s" % (pth, np.asarray(v).ravel()))

print("\n== 5. Pulley (Bifemsh_20mm_Result_pulley_20260926_1652.mat) ==")
d5 = load("Bifemsh_20mm_Result_pulley_20260926_1652.mat")
pb = d5["predBest"]
xbp = np.asarray(d5["xBest"], dtype=float).ravel()
lbp = np.asarray(d5["ctx"].lb, dtype=float).ravel()
ubp = np.asarray(d5["ctx"].ub, dtype=float).ravel()
n = min(xbp.size, lbp.size)
viol = int((xbp[:n] < lbp[:n] - 1e-12).sum() + (xbp[:n] > ubp[:n] + 1e-12).sum())
tzp = np.asarray(pb.TorqueZ, dtype=float).ravel()
tzi = np.asarray(pb.TorqueInsZ, dtype=float).ravel()
oa = np.asarray(pb.offAxisTorque, dtype=float).ravel()
print("fBest=%.5f exitflagP=%s boundViolations=%d infeasibleCount=%s pulleyGain=%s tackleParts=%s optimizeGain=%s" %
      (float(d5["fBest"]), d5["exitflagP"], viol, pb.infeasibleCount, pb.pulleyGain,
       d5["ctx"].tackleLineParts, d5["ctx"].optimizePulleyGain))
print("TorqueZ min %.4f  TorqueInsZ min %.4f  offAxis max %.4f  humanTorque span %.4f..%.4f" %
      (tzp.min(), tzi.min(), oa.max(), huf.min(), huf.max()))
print("gain 81.23->85.74: delta=%.3f (%.2f%%); factors vs 23.16: %.3f / %.3f; offAxis/peak %.3f%%" %
      (abs(tzp.min())-81.2343, (abs(tzp.min())-81.2343)/81.2343*100,
       81.2343/23.1587, abs(tzp.min())/23.1587, oa.max()/abs(tzp.min())*100))

print("\n== 6. BiPulley bifemsh ==")
d6 = load("BiPulley_opensim_bifemsh_r_20260926_1702.mat")
ta = np.asarray(d6["st"].tauAbs, dtype=float)
print("tauTarget=%s fmax=%s (src %s) needPair=%s minPair=%.6f f0=%s fBest=%s exitP=%s" %
      (d6["tauTarget"], d6["spec"].fmax, d6["spec"].fmaxSource, d6["needPair"],
       float(d6["minPair"]), d6["f0"], d6["fBest"], d6["exitflagP"]))
print("st.tauAbs n=%d nan=%d finite=%d maxFinite=%.3f" % (ta.size, int(np.isnan(ta).sum()),
      int((~np.isnan(ta)).sum()), np.nanmax(ta)))
print("spec.bundleMargin =", d6["spec"].bundleMargin, " footprintWidth =", d6["spec"].footprintWidth)
dd = load("BiPulley_opensim_med_gas_r_20260926_1655.mat")
print("med_gas spec.bundleMargin =", dd["spec"].bundleMargin)

print("\n== 7. Gait CSV: means incl. walking/running split ==")
rows = list(csv.DictReader(open(SP + r"\gait_validation_20260930.csv")))
from collections import defaultdict
by = defaultdict(dict)
for r in rows:
    by[r["variant"]][r["reference"]] = float(r["kine_score"])
for v, refs in by.items():
    tune = refs.pop("kine_score")
    vals = np.array(list(refs.values()))
    ong = np.array([s for k, s in refs.items() if k.startswith("ong_")])
    run = np.array([s for k, s in refs.items() if "_Run_" in k])
    worst_ref = min(refs, key=refs.get)
    print("%-11s tune=%.1f | overall=%.1f (n=%d) | ong-only=%.1f (n=%d) | run-only=%.1f (n=%d) | worst=%s %.1f"
          % (v, tune, vals.mean(), vals.size, ong.mean(), ong.size, run.mean(), run.size,
             worst_ref, refs[worst_ref]))

print("\n== 8. Ablation + robustness JSONL ==")
for v in ["s3k", "syn6", "w2lvar"]:
    print("-- prune_results_%s.jsonl WALK rows --" % v)
    for line in open(SP + ("\\prune_results_%s.jsonl" % v)):
        j = json.loads(line)
        if j["mode"] == "WALK":
            kc = j["metrics"].get("kine_components", {})
            extra = ("cyc r/l=%s/%s bilateral=%s knee_min=%s T_r=%s" %
                     (kc.get("n_cycles_r"), kc.get("n_cycles_l"), kc.get("bilateral"),
                      round(kc.get("knee_min", float("nan")), 2), kc.get("T_r"))
                     if kc else "NO-CYCLES")
            print("   %-10s kine=%.4f  %s" % (j["config"], j["metrics"]["kine_score"], extra))
for v in ["s3k", "s3kpruned", "syn6", "w2lvar"]:
    print("-- robust_results_%s.jsonl --" % v)
    fell = []
    for line in open(SP + ("\\robust_results_%s.jsonl" % v)):
        j = json.loads(line)
        m = j["metrics"]
        if j["mode"] == "WALK":
            print("   dty=%-5s WALK kine=%.4f fell=%s tiltMax=%.1f" %
                  (j["dty"], m["kine_score"], m["bal_fell"], m["bal_tilt_max"]))
        elif j["mode"] in ("STAND", "PUSH"):
            if m["bal_fell"]:
                fell.append((j["dty"], j["mode"], j.get("push", "")))
    nsp = sum(1 for line in open(SP + ("\\robust_results_%s.jsonl" % v))
              if json.loads(line)["mode"] in ("STAND", "PUSH"))
    sway = [json.loads(line)["metrics"]["bal_push_sway"]
            for line in open(SP + ("\\robust_results_%s.jsonl" % v))
            if json.loads(line)["mode"] == "PUSH"]
    print("   STAND/PUSH rows n=%d falls=%s maxPushSway=%.4f" % (nsp, fell if fell else "none", max(sway)))

print("\n== 9. Case_40_motion reference ==")
p = SP + r"\gait_refs\Case_40_motion.npz"
if os.path.isfile(p):
    z = np.load(p, allow_pickle=True)
    print("keys:", list(z.keys()))
    for k in z.keys():
        a = np.asarray(z[k])
        if a.size <= 8:
            print(" ", k, "=", a)
        else:
            print(" ", k, "shape", a.shape)
else:
    print("NOT FOUND:", p)
