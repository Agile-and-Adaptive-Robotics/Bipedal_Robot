"""Reconstruct signed muscle knee moments at saved Gait2392 states."""
import argparse
import json
import os
from pathlib import Path
import sys

os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
os.environ['OPENBLAS_NUM_THREADS'] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'functional_gait2392'))

import mujoco
import numpy as np
from mechanics import CoupledArms, JOINTS


def audit(folder, side, end_s, window_s=.12):
    model = mujoco.MjModel.from_binary_path(str(folder / 'model.mjb'))
    state = mujoco.MjData(model)
    arms = CoupledArms(model)
    knee = JOINTS.index('knee_angle_' + side)
    names = [model.actuator(i).name for i in range(model.nu)]
    with np.load(folder / 'aligned_trace.npz') as trace:
        t = trace['t_state']
        selected = np.flatnonzero((t > end_s - window_s) & (t <= end_s))
        torque = np.empty((len(selected), model.nu))
        for row, i in enumerate(selected):
            for dest, source in (('qpos', 'qpos'), ('qvel', 'qvel'),
                                 ('act', 'activation'), ('ctrl', 'control')):
                getattr(state, dest)[:] = trace[source][i]
            mujoco.mj_forward(model, state)
            torque[row] = -arms.update(state)[knee] * state.actuator_force
    if not len(torque):
        raise ValueError('No samples in requested window')
    order = np.argsort(torque.mean(axis=0))
    def rows(ids):
        return [dict(muscle=names[i], mean_Nm=float(torque[:, i].mean()),
                     last_Nm=float(torque[-1, i])) for i in ids]
    return dict(capture=str(folder), side=side, interval_s=[float(t[selected[0]]), float(t[selected[-1]])],
                extension=rows(order[-8:][::-1]), flexion=rows(order[:8]),
                net_muscle_mean_Nm=float(torque.sum(axis=1).mean()),
                net_muscle_last_Nm=float(torque[-1].sum()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('capture', type=Path)
    parser.add_argument('side', choices=('l', 'r'))
    parser.add_argument('end_s', type=float)
    args = parser.parse_args()
    print(json.dumps(audit(args.capture, args.side, args.end_s), indent=2))
