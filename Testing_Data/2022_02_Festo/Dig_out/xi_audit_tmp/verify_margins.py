# verify_margins.py - margins + torque tables from the two optimized mats
import scipy.io as sio
import numpy as np
import os

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"

def margin_report(path, label):
    print(f"\n=== {label}: {os.path.basename(path)} ===")
    m = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
    ks = set(m.keys()) - {'__header__','__version__','__globals__'}
    if 'resultFile' in m:
        print(f" resultFile = {np.asarray(m['resultFile']).flatten()}")
    if 'stamp' in m:
        print(f" stamp = {np.asarray(m['stamp']).flatten()}")
    if 'XiUsed' in m:
        print(f" XiUsed = {np.asarray(m['XiUsed']).flatten()}")
    if 'fBest' in m:
        print(f" fBest = {np.asarray(m['fBest']).flatten()}")
    pb = m.get('predBest')
    tmf = m.get('torqueMarginFraction')
    vht = m.get('validHumanTorque')
    ctx = m.get('ctx')
    if pb is not None and tmf is not None:
        T = np.asarray(pb.Torque).flatten()
        tmf = np.asarray(tmf).flatten()
        phiD = np.asarray(ctx.phiD).flatten() if ctx is not None and hasattr(ctx, 'phiD') else np.arange(1, T.size+1)
        vht = np.asarray(vht).flatten().astype(bool) if vht is not None else np.ones_like(T, bool)
        hu = np.asarray(ctx.humanTorqueAbs).flatten() if ctx is not None and hasattr(ctx, 'humanTorqueAbs') else None
        i = np.nanargmin(np.where(vht, tmf, np.nan))
        print(f" phiD range = [{phiD.min():.3f}, {phiD.max():.3f}] deg, N = {T.size}")
        print(f" MIN torqueMarginFraction = {tmf[i]:.6f} ({tmf[i]*100:.3f}%) at phiD = {phiD[i]:.3f} deg")
        print(f" over VALID frames only: min = {np.nanmin(np.where(vht,tmf,np.nan))*100:.3f}%, mean = {np.nanmean(np.where(vht,tmf,np.nan))*100:.3f}%")
        if hu is not None:
            print(f" robot/human torque at the min-margin frame: {T[i]:.4f} / {hu[i]:.4f} N.m")
            # count violations below target on valid frames
            viol = vht & (T < hu)
            print(f" valid frames with robot < human: {int(viol.sum())} of {int(vht.sum())}")
        # sample table every 10 deg of flexion
        print(" phiD   robot   human   margin%")
        for target in [10, 0, -10, -20, -30, -40, -50, -60, -70, -80, -90, -100, -110, -120]:
            j = int(np.argmin(np.abs(phiD - target)))
            mg = f"{tmf[j]*100:7.2f}" if vht[j] else "     --"
            hu_s = f"{hu[j]:7.3f}" if hu is not None else "      --"
            print(f" {phiD[j]:7.2f} {T[j]:7.3f} {hu_s} {mg}")
    else:
        print(" predBest/torqueMarginFraction missing; keys:", sorted(ks)[:40])

margin_report(os.path.join(BASE, "KneeFlx_20mm_optimized.mat"), "FLEXOR optimized")
margin_report(os.path.join(BASE, "KneeExt_20mm_optimized.mat"), "EXTENSOR optimized")

print("\n=== dated result mats in Results dir (flexor+extensor, Sep) ===")
RES = r"D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results"
for f in sorted(os.listdir(RES)):
    if f.startswith(("Bifemsh_20mm_Result", "Vas_Pam_20mm_Result")):
        print(" ", f)
