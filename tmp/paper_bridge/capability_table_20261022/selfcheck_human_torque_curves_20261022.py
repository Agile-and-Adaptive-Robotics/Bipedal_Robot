"""Self-check (2026-10-22): structure + doctrine checks on
human_torque_curves_20261022.csv, and a peak cross-check against the
legacy requirement curves of record (OpenSim_Vasti_Results.txt,
OpenSim_Bifem_Results.txt) that Mesh_Optimization grades against.

Checks:
 1. 27 rows; 54 columns (3 id cols + 25 theta + 25 tau + caveat).
 2. Every theta row equals the 25-point inclusive uniform grid over
    rom_table_deg[primary_dof] (map echo; tolerance 1e-6 deg).
 3. Every tau column parses finite.
 4. Sign-of-physics spot checks (OpenSim native, knee flexion NEGATIVE):
    vasti group (id 11) tau > 0 (extensors, extension = positive dir);
    each hamstring group (ids 4,5,6,9,10) tau < 0 at mid-RoM;
    iliacus/psoas (26,27) tau > 0 (hip flexion positive dir).
 5. Peak |tau| of the vasti group and of bifemlh_r vs the legacy .sto
    requirement files (informational cross-check, different grid/RoM
    provenance is expected; print both).
Read-only.
"""
import csv
import json
import os

REPO = r'D:\GitHub\Bipedal_Robot'
OUT = os.path.join(REPO, 'tmp', 'paper_bridge', 'capability_table_20261022')
CSV_P = os.path.join(OUT, 'human_torque_curves_20261022.csv')
MAP = os.path.join(REPO, 'Code', 'Matlab', 'Mesh_Optimization',
                   'Routing27_v2', 'actuator_map_27.json')
VASTI = os.path.join(REPO, 'Testing_Data', '2022_02_Festo',
                     'OpenSim_Vasti_Results.txt')
BIFEM = os.path.join(REPO, 'Testing_Data', '2022_02_Festo',
                     'OpenSim_Bifem_Results.txt')
N = 25

with open(MAP) as f:
    cmap = json.load(f)
rom = cmap['rom_table_deg']

rows = []
with open(CSV_P, newline='') as f:
    rdr = csv.reader((ln for ln in f if not ln.startswith('#')))
    header = next(rdr)
    rows = list(rdr)

assert len(header) == 3 + N + N + 1, f'header cols {len(header)}'
assert len(rows) == 27, f'rows {len(rows)}'
print(f'check1 PASS: 27 rows x {len(header)} cols '
      f'(3 + {N} theta + {N} tau + 1 caveat)')

bad_grid = 0
for r in rows:
    aid, _, dof = r[0], r[1], r[2]
    lo, hi = (float(x) for x in rom[dof])
    want = [lo + (hi - lo) * i / (N - 1) for i in range(N)]
    got = [float(x) for x in r[3:3 + N]]
    if any(abs(g - w) > 1e-6 for g, w in zip(got, want)):
        bad_grid += 1
        print(f'  GRID MISMATCH actuator {aid} {dof}')
    taus = [float(x) for x in r[3 + N:3 + 2 * N]]
    assert all(t == t and abs(t) < 1e9 for t in taus), f'nonfinite {aid}'
assert bad_grid == 0
print('check2+3 PASS: all theta grids = 25-pt inclusive uniform over '
      'rom_table_deg; all tau finite')


def tau_of(aid):
    r = next(x for x in rows if x[0] == str(aid))
    return ([float(v) for v in r[3:3 + N]],
            [float(v) for v in r[3 + N:3 + 2 * N]])


th_v, ta_v = tau_of(11)
assert all(t > 0 for t in ta_v), 'vasti should be positive (extension)'
mid = N // 2
for aid in (4, 5, 6, 9, 10):
    _, ta = tau_of(aid)
    assert ta[mid] < 0, f'hamstring actuator {aid} should flex (negative)'
for aid in (26, 27):
    _, ta = tau_of(aid)
    assert ta[mid] > 0, f'hip flexor actuator {aid} should be positive'
print('check4 PASS: sign-of-physics (vasti + extension, hamstrings -, '
      'iliacus/psoas + on hip flexion)')


def sto_peaks(path):
    with open(path, errors='replace') as f:
        lines = f.read().splitlines()
    hdr = next(i for i, ln in enumerate(lines) if ln.strip() == 'endheader')
    cols = [c for c in lines[hdr + 1].replace('\t', ' ').split(' ') if c]
    best = {c: 0.0 for c in cols[2:]}
    for ln in lines[hdr + 2:]:
        if not ln.strip():
            continue
        vals = [float(v) for v in ln.split()]
        for j, c in enumerate(cols[2:], start=2):
            best[c] = max(best[c], abs(vals[j]))
    return cols, best


vc, vp = sto_peaks(VASTI)
bc, bp = sto_peaks(BIFEM)
v_keys = [c for c in vc[2:] if 'vas' in c.lower()]
b_keys = [c for c in bc[2:] if 'bifem' in c.lower() or 'semimem' in c.lower()
          or 'semiten' in c.lower()]
v_group_peak = max(abs(t) for t in ta_v)
_, ta_bf = tau_of(4)
print(f'check5 (informational cross-check vs legacy requirement of record):')
print(f'  vasti group peak here: {v_group_peak:.2f} N*m  | legacy '
      f'per-muscle peaks {[(k, round(vp[k], 2)) for k in v_keys]}'
      f'  legacy group sum peak ~{sum(vp[k] for k in v_keys):.2f}')
print(f'  bifemlh peak here: {max(abs(t) for t in ta_bf):.2f} N*m  | legacy '
      f'{[(k, round(bp[k], 2)) for k in b_keys]}')
print('SELFCHECK DONE')
