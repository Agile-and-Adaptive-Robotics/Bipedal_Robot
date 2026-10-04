"""Probe 2: which model produced OpenSim_Vasti_Results.txt (Opt_run_Ext target)?

Same design as probe_torque_provenance, for vas_int_r/vas_lat_r/vas_med_r
about knee_angle_r. Compares stock gait2392 vs ConnorBipedal_Vastus_Adjusted
(robot-moved vasti routes) against the legacy per-muscle columns.

Run on easteregg2:  D:\\Anaconda\\envs\\opensim\\python.exe probe_vasti_provenance.py
"""
import csv
import opensim as osim

MODELS = {
    'stock': r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_simbody.osim',
    'adjusted': r'D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\ConnorBipedal_Vastus_Adjusted.osim',
}
TARGET = r'D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo\OpenSim_Vasti_Results.txt'
OUT_CSV = r'D:\GitHub\Bipedal_Robot\tmp\scout_probe_20261003\vasti_torque_provenance.csv'
MUSCLES = ['vas_int_r', 'vas_lat_r', 'vas_med_r']
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
        if len(parts) >= 5:
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    return ([r[1] for r in rows], {m: [r[2 + k] for r in rows]
                                   for k, m in enumerate(MUSCLES)})


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
            out[m].append(mu.computeMomentArm(state, coord) * mu.getMaxIsometricForce())
    return angles, out


def interp(angles, taus, grid):
    res = []
    for a in grid:
        if a <= angles[0]:
            res.append(taus[0])
        elif a >= angles[-1]:
            res.append(taus[-1])
        else:
            t = (a - angles[0]) / (angles[-1] - angles[0]) * (N - 1)
            lo = int(t)
            fr = t - lo
            res.append(taus[lo] * (1 - fr) + taus[min(lo + 1, N - 1)] * fr)
    return res


def main():
    leg_ang, leg = legacy_target()
    print(f'legacy vasti: {len(leg_ang)} rows, knee {leg_ang[0]:.1f}..{leg_ang[-1]:.1f} deg')
    for m in MUSCLES:
        print(f'  legacy {m}: {min(leg[m]):.2f}..{max(leg[m]):.2f} N*m')
    results = {}
    for label, path in MODELS.items():
        angles, out = curves(path)
        results[label] = (angles, out)
        for m in MUSCLES:
            print(f'{label:9s} {m}: {min(out[m]):.2f}..{max(out[m]):.2f} N*m')
    print()
    with open(OUT_CSV, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['knee_deg'] + [f'{lb}_{m}' for lb in MODELS for m in MUSCLES]
                   + [f'legacy_{m}' for m in MUSCLES])
        for i in range(N):
            w.writerow([results['stock'][0][i]]
                       + [results[lb][1][m][i] for lb in MODELS for m in MUSCLES]
                       + [leg[m][i] if i < len(leg_ang) else '' for m in MUSCLES])
    for m in MUSCLES:
        rels = {}
        for label in MODELS:
            angles, out = results[label]
            model_i = interp(angles, out[m], leg_ang)
            n = len(model_i)
            mean = sum(abs(x) for x in leg[m]) / n
            rms = (sum((x - y) ** 2 for x, y in zip(model_i, leg[m])) / n) ** 0.5
            rels[label] = rms / mean
            print(f'{label:9s} vs legacy {m}: rms {rms:8.3f} N*m (rel {rms/mean:6.3f})')
        winner = min(rels, key=rels.get)
        print(f'  -> {m}: legacy target matches {winner.upper()}\n')
    print(f'wrote {OUT_CSV}')


if __name__ == '__main__':
    main()
