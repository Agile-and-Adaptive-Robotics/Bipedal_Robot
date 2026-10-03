import scipy.io as sio
import numpy as np

R = r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

def load(f):
    return sio.loadmat(R + "\\" + f, squeeze_me=True, struct_as_record=False)

print("---- Vas_Pam_20mm_Result.mat (extensor design of record) ----")
d = load("Vas_Pam_20mm_Result.mat")
tm = np.asarray(d["torqueMargin"], dtype=float).ravel()
print("torqueMargin: min=%.6f at index %d (0-based) / %d (1-based); n=%d" % (tm.min(), tm.argmin(), tm.argmin()+1, tm.size))
print("xBest full:", np.asarray(d["xBest"]).ravel())
for cand in ["torqueBest", "tauBest", "predBest", "Torque"]:
    if cand in d:
        v = d[cand]
        try:
            print(cand, "->", type(v).__name__, getattr(v, "_fieldnames", None))
        except Exception as e:
            print(cand, e)
pb = d["predBest"]
for f in ["Torque", "TorqueX", "TorqueY", "TorqueZ", "TorqueIns"]:
    if hasattr(pb, f):
        a = np.asarray(getattr(pb, f), dtype=float)
        a = np.nan_to_num(a, nan=np.nan)
        print(f"predBest.{f}: shape {a.shape} min={np.nanmin(a):.4f} max={np.nanmax(a):.4f}")
ctx = d["ctx"]
print("ctx.humanTorque:", np.nanmin(np.asarray(ctx.humanTorque, dtype=float)), np.nanmax(np.asarray(ctx.humanTorque, dtype=float)))
print("ctx.humanTorqueAbs:", np.nanmin(np.asarray(ctx.humanTorqueAbs, dtype=float)), np.nanmax(np.asarray(ctx.humanTorqueAbs, dtype=float)))
print("ctx.humanAngleD:", np.asarray(ctx.humanAngleD).ravel()[:3], "...", np.asarray(ctx.humanAngleD).ravel()[-3:], "n=", np.asarray(ctx.humanAngleD).size)
print("ctx.KMAX:", ctx.KMAX, " ctx.targetPressure:", getattr(ctx, "targetPressure", None), " ctx.BPAcount:", ctx.BPAcount, " ctx.Dia:", getattr(ctx, "Dia", None))
# torque vs human arrays for span check
torq = np.asarray(pb.Torque, dtype=float)
hu = np.asarray(ctx.humanTorque, dtype=float)
n = min(torq.size, hu.size)
print("predBest.Torque span: %.3f .. %.3f  | humanTorque span: %.3f .. %.3f" % (torq.min(), torq.max(), hu.min(), hu.max()))
ang = np.asarray(ctx.humanAngleD, dtype=float).ravel()
print("angle grid: n=%d, from %.3f to %.3f deg, spacing %.4f" % (ang.size, ang.min(), ang.max(), ang[1]-ang[0]))
# release schedule
ri = d["routeInfo"]
print("routeInfo.eliminatedAngleD:", np.asarray(ri.eliminatedAngleD).ravel())
print("routeInfo.eliminationOrder:", np.asarray(ri.eliminationOrder).ravel())
print("routeInfo.activationOrder:", np.asarray(ri.activationOrder).ravel())
print("routeInfo.activeAtFullFlexion:", np.asarray(ri.activeAtFullFlexion).ravel())

print()
print("---- Vas_Pam_20mm_Result_20260925.mat ----")
d2 = load("Vas_Pam_20mm_Result_20260925.mat")
tm2 = np.asarray(d2["torqueMargin"], dtype=float).ravel()
print("torqueMargin min %.6f at idx %d (0-based)" % (tm2.min(), tm2.argmin()))
print("xBest full:", np.asarray(d2["xBest"]).ravel())
ri2 = d2["routeInfo"]
print("eliminatedAngleD:", np.asarray(ri2.eliminatedAngleD).ravel())
pb2 = d2["predBest"]
t2 = np.asarray(pb2.Torque, dtype=float)
hu2 = np.asarray(d2["ctx"].humanTorque, dtype=float)
print("Torque span %.3f..%.3f | humanTorque span %.3f..%.3f" % (t2.min(), t2.max(), hu2.min(), hu2.max()))

print()
print("---- Bifemsh_20mm_Result.mat (flexor) ----")
d3 = load("Bifemsh_20mm_Result.mat")
print("fBest:", d3["fBest"], " Jseed:", d3["Jseed"])
tmf = np.asarray(d3["torqueMarginFraction"], dtype=float).ravel()
print("torqueMarginFraction min %.6f at idx %d; angles at idx: humanAngleD[%d]=" % (tmf.min(), tmf.argmin(), tmf.argmin()), end="")
print(np.asarray(d3["ctx"].humanAngleD, dtype=float).ravel()[tmf.argmin()])
pb3 = d3["predBest"]
for f in ["Torque", "TorqueX", "TorqueY", "TorqueZ", "offAxisTorque", "TorqueIns"]:
    if hasattr(pb3, f):
        a = np.asarray(getattr(pb3, f), dtype=float)
        if a.size:
            print(f"predBest.{f}: shape {a.shape} min={np.nanmin(a):.4f} max={np.nanmax(a):.4f}")
hu3 = np.asarray(d3["ctx"].humanTorque, dtype=float)
print("humanTorque span: %.4f .. %.4f" % (hu3.min(), hu3.max()))
print("ctx.KMAX:", d3["ctx"].KMAX, "BPAcount:", d3["ctx"].BPAcount, "Dia:", getattr(d3["ctx"], "Dia", None), "targetPressure:", getattr(d3["ctx"], "targetPressure", None))

print()
print("---- Bifemsh_20mm_Result_pulley_20260926_1652.mat ----")
d4 = load("Bifemsh_20mm_Result_pulley_20260926_1652.mat")
print("fBest:", d4["fBest"], " exitflagP:", d4["exitflagP"], " Jseed:", d4.get("Jseed"))
pb4 = d4["predBest"]
for f in ["Torque", "TorqueX", "TorqueY", "TorqueZ", "offAxisTorque", "TorqueIns", "TorqueInsZ", "pulleyGain"]:
    if hasattr(pb4, f):
        a = np.asarray(getattr(pb4, f), dtype=float)
        if a.size:
            print(f"predBest.{f}: shape {a.shape} min={np.nanmin(a):.4f} max={np.nanmax(a):.4f}")
hu4 = np.asarray(d4["ctx"].humanTorque, dtype=float)
print("humanTorque span: %.4f .. %.4f" % (hu4.min(), hu4.max()))

print()
print("---- BiPulley_opensim_bifemsh_r_20260926_1702.mat ----")
d5 = load("BiPulley_opensim_bifemsh_r_20260926_1702.mat")
print("keys:", [k for k in d5 if not k.startswith("__")])
st = d5["st"]
print("st fields:", st._fieldnames)
ta = np.asarray(st.tauAbs, dtype=float)
print("st.tauAbs shape:", ta.shape, "n=", ta.size, "allNaN=", bool(np.isnan(ta).all()), "nanCount=", int(np.isnan(ta).sum()))
print("st.obj:", st.obj, " st.worstC:", st.worstC, " st.margin:", st.margin)
print("fBest:", d5["fBest"], "exitflagP:", d5["exitflagP"], "exitflagS:", d5["exitflagS"])
# initial objective: look for surrogate/log fields
for k in d5:
    if k.startswith("__"): continue
    v = d5[k]
    if isinstance(v, np.ndarray) and v.dtype == object:
        pass
spec = d5.get("spec", None)
if spec is not None:
    print("spec fields:", spec._fieldnames if hasattr(spec, "_fieldnames") else type(spec))
    for f in ["targetTorque", "MIF", "targetName", "muscle"]:
        if hasattr(spec, f): print("spec.%s =" % f, getattr(spec, f))
for k in [k for k in d5 if not k.startswith("__")]:
    if k in ("st","spec","ctx","routeCtx"): continue
    v = np.asarray(d5[k])
    if v.size <= 10:
        print(k, "=", v.ravel())
