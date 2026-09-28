# confirmer_dump_20260926.py - independent re-derivation of every quoted mat number
# (Xi/GoF confirmer pass). Read-only: no mats written.
import scipy.io as sio
import numpy as np
import os, datetime

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"
RES = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"

def mt(p):
    return datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M:%S")

def flat(v):
    return np.asarray(v).flatten()

print("== 1. flexor front: minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat ==")
S = sio.loadmat(os.path.join(BASE, "minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat"),
                squeeze_me=True, struct_as_record=False)
xC = flat(S['xCols']).astype(int)
print(" xCols =", xC)
for r in (77, 107, 1):
    print(f" filtered_results({r}, xCols) =", ["%.6g" % x for x in flat(S['filtered_results'][r-1, xC-1])])

print("\n== 2. extensor front: minimizeExt10mmX3_results_20260920_flxr77.mat ==")
E = sio.loadmat(os.path.join(BASE, "minimizeExt10mmX3_results_20260920_flxr77.mat"),
                squeeze_me=True, struct_as_record=False)
xCe = flat(E['xCols']).astype(int)
print(" xCols =", xCe)
for r in (32,):
    print(f" filtered_results({r}, xCols) =", ["%.6g" % x for x in flat(E['filtered_results'][r-1, xCe-1])])
for k in ('LOCKSRC', 'KMODE', 'BOUNDS', 'flxRow', 'xi1lock', 'xi2lock'):
    if k in E:
        print(f" {k} =", flat(E[k]))
if 'baselineScores' in E:
    bs = np.asarray(E['baselineScores'])
    print(" baselineScores shape =", bs.shape)
    print(" baselineScores =", np.array2string(bs, precision=4, max_line_width=200))

print("\n== 3. old extensor front: minimizeExt10mmX3_results_20260910_noT3.mat row 1 ==")
E0 = sio.loadmat(os.path.join(BASE, "minimizeExt10mmX3_results_20260910_noT3.mat"),
                 squeeze_me=True, struct_as_record=False)
xC0 = flat(E0['xCols']).astype(int)
print(" xCols =", xC0, " row 1 =", ["%.6g" % x for x in flat(E0['filtered_results'][0, xC0-1])])

print("\n== 4. Bifemsh_20mm_Result.mat (flexor design of record claim) ==")
p = os.path.join(RES, "Bifemsh_20mm_Result.mat")
print(" mtime =", mt(p))
B = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
for k in ('resultFile', 'stamp', 'XiUsed', 'fBest'):
    if k in B:
        print(f" {k} =", flat(B[k]))
tmf = B.get('torqueMarginFraction'); pb = B.get('predBest'); ctx = B.get('ctx')
if tmf is not None and pb is not None:
    t = flat(tmf); T = flat(pb.Torque)
    phiD = flat(ctx.phiD) if ctx is not None and hasattr(ctx, 'phiD') else np.arange(1, t.size+1)
    vht = flat(B['validHumanTorque']).astype(bool) if 'validHumanTorque' in B else np.ones(t.size, bool)
    i = int(np.nanargmin(np.where(vht, t, np.nan)))
    print(f" min torqueMarginFraction = {t[i]:.6f} ({t[i]*100:.3f}%) at phiD = {phiD[i]:.3f} deg (valid-only argmin)")
    j = int(np.nanargmin(t))
    print(f" min over ALL frames      = {t[j]:.6f} ({t[j]*100:.3f}%) at phiD = {phiD[j]:.3f} deg")

print("\n== 5. Vas_Pam_20mm_Result.mat (extensor base, claimed overwritten 09-25) ==")
p = os.path.join(RES, "Vas_Pam_20mm_Result.mat")
print(" mtime =", mt(p))
V = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
ks = sorted(k for k in V.keys() if not k.startswith('__'))
print(" keys:", ks[:40])
for k in ('liveRun', 'fBest', 'minMargin', 'stamp', 'XiUsed', 'resultFile'):
    if k in V:
        print(f" {k} =", flat(V[k]))

print("\n== 6. Vas_Pam_20mm_Result_20260925.mat (09-25 re-run) ==")
p = os.path.join(RES, "Vas_Pam_20mm_Result_20260925.mat")
print(" mtime =", mt(p))
V2 = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
ks2 = sorted(k for k in V2.keys() if not k.startswith('__'))
print(" keys:", ks2[:40])
for k in ('liveRun', 'fBest', 'minMargin', 'removedText', 'transitionIdx'):
    if k in V2:
        try:
            s = flat(V2[k])
            print(f" {k} =", s if s.dtype.kind in 'US' else np.array2string(s, precision=5))
        except Exception as ex:
            print(f" {k}: <{ex}>")
if 'routeInfo' in V2:
    ri = V2['routeInfo']
    for fn in ri._fieldnames:
        try:
            v = flat(getattr(ri, fn))
            if v.size <= 12:
                print(f" routeInfo.{fn} =", np.array2string(v, precision=4, max_line_width=200))
        except Exception:
            pass

print("\n== 7. static tables: FlxPinBPASet / ExtPinBPASet ==")
m = sio.loadmat(os.path.join(BASE, "FlxPinBPASet.mat"), squeeze_me=True, struct_as_record=False)
kf = np.atleast_1d(m['kf'])
for i, s in enumerate(kf, 1):
    ak = flat(s.Ak)
    print(f" kf({i}): rest={float(s.rest):.4f} ten={float(s.ten):.4f} P={float(s.P):.1f} Fm={float(s.Fm):.1f} "
          f"dBPA={float(s.dBPA)*1000:.1f}mm Ak=[{ak.min():.1f},{ak.max():.1f}] M(Loc rows)={flat(s.Loc).shape[0]}")
e = sio.loadmat(os.path.join(BASE, "ExtPinBPASet.mat"), squeeze_me=True, struct_as_record=False)
ke = np.atleast_1d(e['ke'])
for i, s in enumerate(ke, 1):
    ak = flat(s.Ak)
    print(f" ke({i}): rest={float(s.rest):.4f} ten={float(s.ten):.4f} P={float(s.P):.1f} Fm={float(s.Fm):.1f} "
          f"dBPA={float(s.dBPA)*1000:.1f}mm Ak=[{ak.min():.1f},{ak.max():.1f}] M(Loc rows)={flat(s.Loc).shape[0]}")

print("\nCONFIRMER DUMP DONE (read-only)")
