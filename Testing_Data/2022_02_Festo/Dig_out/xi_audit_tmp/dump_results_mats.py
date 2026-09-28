# dump_results_mats.py - XiUsed/fBest/margins from base + latest dated result mats
import scipy.io as sio
import numpy as np
import os, datetime

RES = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
files = [
    "Bifemsh_20mm_Result.mat",
    "Bifemsh_20mm_Result_20260920_1443.mat",
    "Bifemsh_20mm_Result_20260921_2037.mat",
    "Vas_Pam_20mm_Result.mat",
    "Vas_Pam_20mm_Result_20260920_1519.mat",
    "Vas_Pam_20mm_Result_20260925.mat",
]
for f in files:
    p = os.path.join(RES, f)
    if not os.path.isfile(p):
        print(f"{f}: MISSING"); continue
    mt = datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M")
    try:
        m = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
    except Exception as ex:
        print(f"{f} (mtime {mt}): loadmat FAIL {type(ex).__name__}: {ex}"); continue
    print(f"\n--- {f} (mtime {mt}) ---")
    for k in ['XiUsed', 'fBest', 'stamp', 'resultFile']:
        if k in m:
            print(f"  {k} = {np.asarray(m[k]).flatten()}")
    pb = m.get('predBest'); tmf = m.get('torqueMarginFraction'); ctx = m.get('ctx')
    if pb is not None and tmf is not None and ctx is not None:
        T = np.asarray(pb.Torque).flatten()
        tmf = np.asarray(tmf).flatten()
        phiD = np.asarray(ctx.phiD).flatten()
        vht = np.asarray(m['validHumanTorque']).flatten().astype(bool) if 'validHumanTorque' in m else np.ones(T.size, bool)
        hu = np.asarray(ctx.humanTorqueAbs).flatten() if hasattr(ctx, 'humanTorqueAbs') else None
        if vht.any():
            i = int(np.nanargmin(np.where(vht, tmf, np.nan)))
            print(f"  MIN margin = {np.where(vht,tmf,np.nan)[i]*100:.3f}% at phiD {phiD[i]:.3f} deg"
                  + (f" (robot {T[i]:.3f} vs human {hu[i]:.3f} N.m)" if hu is not None else ""))
        print(f"  torqueMarginFraction: min {np.nanmin(np.where(vht,tmf,np.nan))*100:.3f}%  mean(valid) {np.nanmean(np.where(vht,tmf,np.nan))*100:.3f}%")
    elif 'Torque' in m or 'Torque0' in m or 'Bifemsh_T' in m:
        for k in ['Angle', 'Torque', 'Torque0', 'Torque2', 'Torque3', 'Bifemsh_T', 'Bifemsh_MA', 'H', 'Xi0', 'Xi1', 'Xi2', 'Xi3']:
            if k in m:
                a = np.asarray(m[k])
                print(f"  {k}: shape {a.shape}" + (f" = {a.flatten()}" if a.size <= 4 else ""))
