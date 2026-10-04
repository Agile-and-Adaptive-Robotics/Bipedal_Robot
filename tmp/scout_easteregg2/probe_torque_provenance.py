"""Probe: which model produced OpenSim_Bifem_Results.txt (the Opt_run human target)?

Computes max-isometric knee torque curves (tau = MIF * moment arm) for
bifemlh_r / bifemsh_r from (a) stock gait2392_simbody.osim and (b)
ConnorBipedal_Bifemsh_Adjusted.osim, then compares both against the legacy
target file column-by-column. Read-only; writes one CSV + prints a verdict.

Run on easteregg2:  D:\\Anaconda\\envs\\opensim\\python.exe probe_torque_provenance.py
"""
import csv
import opensim as osim

MODELS = {
    'stock': r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_simbody.osim',
    'adjusted': r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\ConnorBipedal_Bifemsh_Adjusted.osim',
}
TARGET = r'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo\OpenSim_Bifem_Results.txt'
OUT_CSV = r'D:\GitHub\Bipedal_Robot\tmp\scout_probe_20261003\bifem_torque_provenance.csv'
MUSCLES = ['bifemlh_r', 'bifemsh_r']
QNAME = 'knee_angle_r'
N = 100
DEG_MIN, DEG_MAX = -120.0, 10.0


def legacy_target():
    rows = []
    with open(TARGET) as f:
        lines = f.read().splitlines()
    hdr = next(i for i, ln in enumerate(lines) if ln.strip() == 'endheader')
    for ln in lines[hdr + 1:]:
        parts = ln.split()
        if len(parts) >= 4:
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue  # column-name row
    # col1 time(idx), col2 knee angle deg, col3 bifemlh_r, col4 bifemsh_r
    return ([r[1] for r in rows], [r[2] for r in rows], [r[3] for r in rows])


def curves(model_path):
    model = osim.Model(model_path)
    state = model.initSystem()
    coord = model.getCoordinateSet().get(QNAME)
    out = {m: [] for m in MUSCLES}
    angles = []
    ms = model.getMuscles()
    for i in range(N):
        deg = DEG_MIN + (DEG_MAX - DEG_MIN) * i / (N - 1)
        coord.setValue(state, deg * 3.141592653589793 / 180.0)
        model.realizePosition(state)
        angles.append(deg)
        for m in MUSCLES:
            mu = ms.get(m)
            arm = mu.computeMomentArm(state, coord)      # meters
            mif = mu.getMaxIsometricForce()
            out[m].append(arm * mif)                     # N*m, flexion-negative
    return angles, out


def cmp(a, b):
    n = min(len(a), len(b))
    sa, sb = a[:n], b[:n]
    mean = sum(abs(x) for x in sb) / n
    rms = (sum((x - y) ** 2 for x, y in zip(sa, sb)) / n) ** 0.5
    return rms, rms / mean if mean else float('inf')


def main():
    leg_ang, leg_lh, leg_sh = legacy_target()
    print(f'legacy target: {len(leg_ang)} rows, knee {leg_ang[0]:.1f}..{leg_ang[-1]:.1f} deg')
    print(f'legacy bifemlh_r range {min(leg_lh):.3f}..{max(leg_lh):.3f}, '
          f'bifemsh_r {min(leg_sh):.3f}..{max(leg_sh):.3f} N*m')
    results = {}
    for label, path in MODELS.items():
        angles, out = curves(path)
        results[label] = (angles, out)
        for m in MUSCLES:
            print(f'{label:9s} {m}: tau range {min(out[m]):.3f}..{max(out[m]):.3f} N*m')
    print()
    verdict = {}
    for label in MODELS:
        angles, out = results[label]
        for m, leg in (('bifemlh_r', leg_lh), ('bifemsh_r', leg_sh)):
            # interpolate model curve onto the legacy angle grid
            taus = out[m]
            interp = []
            for a in leg_ang:
                if a <= angles[0]:
                    interp.append(taus[0])
                elif a >= angles[-1]:
                    interp.append(taus[-1])
                else:
                    t = (a - angles[0]) / (angles[-1] - angles[0]) * (N - 1)
                    lo = int(t)
                    frac = t - lo
                    interp.append(taus[lo] * (1 - frac) + taus[min(lo + 1, N - 1)] * frac)
            rms, rel = cmp(interp, leg)
            verdict[(label, m)] = (rms, rel)
            print(f'{label:9s} vs legacy {m}: rms {rms:8.3f} N*m  (rel {rel:6.3f})')
    print()
    for m in MUSCLES:
        s = verdict[('stock', m)]
        a = verdict[('adjusted', m)]
        winner = 'STOCK' if s < a else 'ADJUSTED'
        print(f'{m}: legacy target matches {winner} '
              f'(stock rel {s[1]:.4f} vs adjusted rel {a[1]:.4f})')
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['knee_deg'] + [f'{lb}_{m}_tau_Nm' for lb in MODELS for m in MUSCLES]
                   + ['legacy_bifemlh_r', 'legacy_bifemsh_r'])
        for i in range(N):
            w.writerow([results['stock'][0][i]]
                       + [results[lb][1][m][i] for lb in MODELS for m in MUSCLES]
                       + ['', ''])
        for i, a in enumerate(leg_ang):
            w.writerow([a, '', '', '', '', leg_lh[i], leg_sh[i]])
    print(f'wrote {OUT_CSV}')


if __name__ == '__main__':
    main()
