"""Resolve the grac_r step-magnitude discrepancy between
_probe_osim_limits_2d_20261023.py (7.0573% of column peak) and
_probe_step_locations_20261023.py (2.04%): compute the knee-angle column of
moment arms at several hip flexion values THREE ways per column:
  A) computeMomentArm only (location-probe sequence),
  B) computeMomentArm + getLength + isActive interleaved (main-probe seq),
  C) fresh initSystem per column, sequence A.
Print per-column max |d(arm)|, column peak, and ratio. Pinned 4.4.1.
"""
import math
import os
import sys

import opensim as osim

MODEL = os.path.join(r'D:\GitHub\Bipedal_Robot', 'Solid_Models', 'OpenSim',
                     'Gait2392_Robotbody', 'gait2392_simbody.osim')
DEG2RAD = math.pi / 180.0
HIPS = [-15.0, 0.0, 11.0, 45.0, 80.0, 110.0]


def col(model, muscle, hip_deg, interleave, fresh):
    if fresh:
        state = model.initSystem()
    else:
        state = model.states()  # placeholder, replaced below
    return None


def main():
    if str(osim.__version__) != '4.4.1':
        sys.exit('PIN VIOLATION')
    model = osim.Model(MODEL)
    base = model.initSystem()
    coords = model.getCoordinateSet()
    mu = model.getMuscles().get('grac_r')
    knee = coords.get('knee_angle_r')
    hip = coords.get('hip_flexion_r')
    k_degs = [-120.0 + i for i in range(131)]

    def run(hip_deg, interleave, state):
        hip.setValue(state, hip_deg * DEG2RAD)
        arms = []
        for k in k_degs:
            knee.setValue(state, k * DEG2RAD)
            model.realizePosition(state)
            arms.append(mu.computeMomentArm(state, knee))
            if interleave:
                mu.getLength(state)
                pts = mu.getGeometryPath().getPathPointSet()
                for i in range(pts.getSize()):
                    pts.get(i).isActive(state)
        peak = max(abs(v) for v in arms)
        dmax = max(abs(arms[i + 1] - arms[i]) for i in range(len(arms) - 1))
        return peak, dmax

    print('hip    A(arm-only,shared) B(interleave,shared) '
          'C(arm-only,fresh-per-col)')
    for h in HIPS:
        st = model.initSystem()
        pa, da = run(h, False, st)
        st = model.initSystem()
        pb, db = run(h, True, st)
        st = model.initSystem()
        pc, dc = run(h, False, st)
        print(f'{h:6.1f}  {da/pa*100:8.4f}%       {db/pb*100:8.4f}%'
              f'         {dc/pc*100:8.4f}%   '
              f'(peaks {pa:.5f} {pb:.5f} {pc:.5f}; d {da:.6f} {db:.6f} '
              f'{dc:.6f})')

    # determinism: repeat column hip=11 twice, same sequence
    st = model.initSystem()
    p1, d1 = run(11.0, False, st)
    st = model.initSystem()
    p2, d2 = run(11.0, False, st)
    print(f'repeat hip=11: {d1/p1*100:.6f}% vs {d2/p2*100:.6f}% '
          f'identical={abs(d1/p1 - d2/p2) < 1e-12}')
    print('DONE')


if __name__ == '__main__':
    main()
