"""Repeat the original Gait2392 static activation calculation with four compliant tendons.

Supply bsolve_ik.py --fk-signs output. Its complete poses, inverse joint torques,
finite-difference moment arms, zero muscle speeds, and activation-penalty rule
are preserved. The isolated frames have no prescribed fiber-velocity history.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import mujoco
import numpy as np
from scipy.optimize import lsq_linear, minimize

from coupled_muscle_context import SPINAL, FIDELITY
sys.path[:0] = [str(SPINAL), str(FIDELITY)]
import bsolve_ik as original
from selective_force import MODEL, MUSCLES
from thelen_dynamic import DynamicThelenBank, parameters_from_osim


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    path, out = args.input.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    after = out / 'bsolve_out.npz'
    if after.exists():
        parser.error('An activation solution already exists in the output directory')
    with np.load(path, allow_pickle=True) as frozen:
        recording = {key: frozen[key] for key in frozen.files}
    model, data = original.build_air_model(id_mode=True)
    names = recording['act_names'].tolist()
    assert names == original.act_names_of(model)
    joints = recording['fit_joints'].tolist()
    assert joints == list(original.FIT_JOINTS)
    rows = np.array([model.joint(name).dofadr[0] for name in joints])
    arms = -original.fd_moments(model, data, recording['qpos'], rows)
    selected = np.array([names.index(n) for n in MUSCLES])
    other = np.setdiff1d(np.arange(model.nu), selected)
    bank = DynamicThelenBank(parameters_from_osim(MODEL, MUSCLES))
    before = recording['acts'].copy()
    force_trace = np.empty_like(before)
    identity = np.eye(len(other))
    max_native_difference = 0.
    max_equilibrium_error = 0.
    max_optimality = 0.
    alternate_solver_frames = []

    for k, position in enumerate(recording['qpos']):
        data.qpos[:] = position
        data.qvel[:] = 0.
        data.act[:] = 0.
        mujoco.mj_forward(model, data)
        passive = -data.actuator_force.copy()
        np.testing.assert_allclose(data.actuator_length, recording['L_mj'][k], atol=1e-12, rtol=0)
        lengths = data.actuator_length[selected].copy()
        data.act[:] = 1.
        mujoco.mj_forward(model, data)
        capacity = -data.actuator_force - passive
        rho = arms[k]
        target = recording['tau_fit'][k]

        # Rebuild the original optimization as a per-frame positive control for sign, moment-arm, and force conventions.
        matrix = rho * capacity
        norm = np.linalg.norm(matrix, axis=0)
        good = norm > 1e-9
        regularization = .05 * np.median(norm[good])
        native = lsq_linear(np.vstack([matrix, regularization * np.eye(model.nu)]),
                            np.r_[target - rho @ passive, np.zeros(model.nu)],
                            bounds=(0., 1.), tol=1e-10, max_iter=200)
        verification = np.where(good, native.x, 0.)
        max_native_difference = max(max_native_difference, np.max(np.abs(verification - before[k])))

        def tension(command):
            return bank.evaluate(lengths, command, dynamic=False).copy()

        passive[selected] = tension(np.zeros(4))
        capacity[selected] = tension(np.ones(4)) - passive[selected]
        norm = np.linalg.norm(rho * capacity, axis=0)
        regularization = .05 * np.median(norm[norm > 1e-9])
        matrix = np.vstack([rho[:, other] * capacity[other], regularization * identity])
        native_passive = rho[:, other] @ passive[other]

        def remaining(command):
            tension_four = tension(command)
            rhs = target - native_passive - rho[:, selected] @ tension_four
            sol = lsq_linear(matrix, np.r_[rhs, np.zeros(len(other))],
                             bounds=(0., 1.), method='bvls', tol=1e-10, max_iter=200)
            if not sol.success:
                raise RuntimeError(f'Native actuator solve failed at {k}')
            activation = np.empty(model.nu)
            if np.any(sol.x < -1e-12) or np.any(sol.x > 1. + 1e-12):
                raise RuntimeError(f'Native solver returned an infeasible activation at {k}')
            # BVLS can represent an active bound one ulp beyond [0, 1].
            activation[other], activation[selected] = np.clip(sol.x, 0., 1.), command
            forces = passive + capacity * activation
            forces[selected] = tension_four
            return activation, forces

        def objective(command):
            activations, forces = remaining(command)
            residual = rho @ forces - target
            upper, lower = np.minimum(1., command + 1e-5), np.maximum(.01, command - 1e-5)
            slope = (tension(upper) - tension(lower)) / (upper - lower)
            gradient = slope * (rho[:, selected].T @ residual) + regularization**2 * command
            return .5 * (residual @ residual + regularization**2 * (activations @ activations)), gradient

        best = minimize(objective, np.clip(before[k, selected], .01, 1.), jac=True, method='L-BFGS-B',
                        bounds=[(.01, 1.)] * 4, options={'maxiter': 160, 'maxls': 40, 'ftol': 1e-12, 'gtol': 1e-7})
        if not best.success:
            alternate_solver_frames.append(k)
            best = minimize(objective, best.x, jac=True, method='SLSQP', bounds=[(.01, 1.)] * 4,
                            options={'maxiter': 200, 'ftol': 1e-10})
            if not best.success:
                raise RuntimeError(f'Four-muscle recruitment failed at frame {k}: {best.message}')
        activation, fitted_force = remaining(best.x)
        if np.any(bank.bounds):
            raise RuntimeError(f'Fiber equilibrium bound at frame {k}')
        # Source muscles have a minimum activation of 0.01, including at a zero command.
        recording['acts'][k] = activation
        force_trace[k] = fitted_force
        recording['res_norm'][k] = np.linalg.norm(rho @ fitted_force - target) / max(np.linalg.norm(target - rho @ passive), 1e-9)
        max_equilibrium_error = max(max_equilibrium_error, float(np.max(np.abs(bank.error))))
        max_optimality = max(max_optimality, float(np.max(np.abs(best.x - np.clip(best.x - best.jac, .01, 1.)))))

    if max_native_difference > 1e-6 or max_equilibrium_error > 1e-5:
        raise RuntimeError(f'Native/in-series equilibrium verification failed: {max_native_difference}, {max_equilibrium_error}')
    t = recording['t']
    t0, t1 = recording['hs_r'][:2]
    keep = (t >= t0) & (t < t1)
    bins = np.clip(((t[keep] - t0) / (t1 - t0) * 20).astype(int), 0, 19)
    for suffix, side in [('', 'r'), ('_l', 'l')]:
        label = 'grp_prof' + suffix
        for j, group in enumerate(recording[label + '_names']):
            indices = [i for i, name in enumerate(names)
                       if (info := original.classify(name)) and info.side == side and group in info.groups and recording['Fmax'][i] > 5]
            recording[label][j] = [recording['acts'][keep][bins == i][:, indices].mean() for i in range(20)]
    np.savez_compressed(after, **recording)
    np.savez_compressed(out / 'fitted_tension.npz', time=t, act_names=names, force_N=force_trace)
    metrics = {'source': str(path), 'native_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'maximum_native_activation_difference': max_native_difference,
               'maximum_fiber_equilibrium_error_N': max_equilibrium_error,
               'maximum_projected_first_order_residual': max_optimality,
               'alternate_solver_frames': alternate_solver_frames,
               'rms_change_in_activation': float(np.sqrt(np.mean((recording['acts'] - before)**2))),
               'mean_inverse_residual': float(recording['res_norm'].mean()),
               'muscles': MUSCLES, 'frame_count': len(t), 'muscle_velocity': 'zero, as in the original static force step',
               'geometry': 'Anatomically signed input poses, same inverse torques, same dependent-joint finite differences',
               'four_muscle_activation_bounds': [.01, 1.], 'remaining_activation_bounds': [0., 1.]}
    (out / 'recruitment.json').write_text(json.dumps(metrics, indent=2), encoding='utf-8')
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
