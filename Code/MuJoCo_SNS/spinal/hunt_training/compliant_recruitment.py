"""Recalculate gait recruitment with bilateral compliant-tendon RF and BFSH.

Run prepare_targets.py first to calculate the target torques, moment arms, and
native capacities. Hold those data fixed. Only four of the 92 force equations
change; activations for all 92 are solved again. Commands are held over each
sample interval, and path lengths are interpolated at the fiber integration step.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from numba import njit
from scipy.optimize import lsq_linear, minimize

SPINAL = Path(__file__).resolve().parents[1]
FIDELITY = SPINAL / 'physics_fidelity'
sys.path.insert(0, str(FIDELITY))
from thelen_dynamic import evaluate_one, parameters_from_osim

MUSCLES = ('rect_fem_r', 'rect_fem_l', 'bifemsh_r', 'bifemsh_l')
DEFAULT_MODEL = SPINAL.parents[2] / 'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim'


@njit(cache=True)
def interval(previous, start, end, activation, parameters, elapsed, max_step):
    steps = max(1, int(np.ceil(elapsed / max_step)))
    dt = elapsed / steps if elapsed > 0 else max_step
    state = previous.copy()
    force = np.zeros(len(activation))
    maximum_error = 0.
    for j in range(steps):
        lengths = start + (end - start) * ((j + 1) / steps)
        for i in range(len(activation)):
            s, f, error, bound, _, _ = evaluate_one(
                state[i], lengths[i], activation[i], parameters[i], dt, elapsed > 0)
            if bound or not np.isfinite(f):
                raise ValueError('Fiber equilibrium outside numerical bounds')
            maximum_error = max(maximum_error, abs(error))
            state[i], force[i] = s, f
    return state, force, maximum_error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True, help='Output directory from prepare_targets.py')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--osim', type=Path, default=DEFAULT_MODEL)
    parser.add_argument('--max-step', type=float, default=.002)
    args = parser.parse_args()
    if args.max_step <= 0 or not np.isfinite(args.max_step):
        parser.error('Fiber integration step must be positive and finite')
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if (out/'targets.npz').exists() or out.resolve() == args.native.resolve():
        parser.error('Use a fresh output directory to preserve the reference')
    source = args.native/'targets.npz'
    with np.load(source, allow_pickle=False) as native:
        z = {k: native[k] for k in native.files}
    meta = json.loads((args.native/'metadata.json').read_text())
    params = parameters_from_osim(args.osim, MUSCLES)
    ids = np.array([z['act_names'].tolist().index(m) for m in MUSCLES])
    nmus = len(z['act_names'])
    prev_state = params[:,1].copy()
    start = z['length'][0,ids].copy()
    other = np.setdiff1d(np.arange(nmus), ids)
    eye = np.eye(len(other))
    max_error = 0.
    optimality = []
    for i, time in enumerate(z['t']):
        elapsed = float(time-z['t'][i-1]) if i else 0.
        end = z['length'][i,ids].copy()
        def evaluate(a):
            return interval(prev_state, start, end, a, params, elapsed, args.max_step)
        rho = z['rho'][i]
        target = z['torque'][i]
        f0 = z['passive_force'][i].copy()
        capacity = z['unit_force'][i].copy()
        _, f0[ids], _ = evaluate(np.zeros(len(ids)))
        _, full, _ = evaluate(np.ones(len(ids)))
        capacity[ids] = full - f0[ids]
        norms = np.linalg.norm(rho*capacity,axis=0)
        lam = .025 * np.median(norms[norms>1e-6])
        linear_matrix = np.vstack((rho[:,other]*capacity[other], lam*eye))
        fixed_torque = rho[:,other]@f0[other]
        def fit_remaining(a4):
            f4 = evaluate(a4)[1]
            rhs = target - fixed_torque - rho[:,ids]@f4
            fit = lsq_linear(linear_matrix, np.r_[rhs, np.zeros(len(other))],
                             bounds=(0.,.98), method='bvls', tol=1e-8, max_iter=160)
            if not fit.success:
                raise RuntimeError(f'Native recruitment failed at {time} s')
            a = np.empty(nmus)
            a[ids], a[other] = a4, fit.x
            f = f0 + capacity*a
            f[ids] = f4
            return a, f
        def objective(a4):
            a, f = fit_remaining(a4)
            err = rho@f-target
            hi, lo = np.minimum(1., a4+1e-5), np.maximum(.01, a4-1e-5)
            slope = (evaluate(hi)[1]-evaluate(lo)[1])/(hi-lo)
            gradient = (rho[:,ids].T@err)*slope + lam**2*a4
            return .5*(err@err+lam**2*(a@a)), gradient
        guess = np.clip(z['activation'][i,ids],.01,.98)
        sol = minimize(objective, guess, method='L-BFGS-B', jac=True, bounds=[(.01,.98)]*len(ids),
                       options={'maxiter':150,'maxls':40,'ftol':1e-11,'gtol':1e-6})
        if not sol.success:
            raise RuntimeError(f'Recruitment failed at {time} s: {sol.message}')
        z['activation'][i], f = fit_remaining(sol.x)
        z['force'][i] = f
        z['passive_force'][i], z['unit_force'][i] = f0, capacity
        z['residual'][i] = np.linalg.norm(rho@f-target)/max(np.linalg.norm(target-rho@f0),1e-6)
        prev_state, f_check, error = evaluate(sol.x)
        np.testing.assert_allclose(f[ids],f_check,rtol=0,atol=0)
        start = end
        max_error = max(max_error,error)
        optimality.append(float(np.max(np.abs(sol.x - np.clip(sol.x-sol.jac,.01,.98)))))
    if max_error > 1e-5:
        raise RuntimeError(f'Fiber balance residual {max_error} N')
    for side_index, side in enumerate(('r','l')):
        for feature in (4,5):
            members = meta['feature_members'][side][feature]
            z['features'][:,side_index,feature] = np.maximum(
                z['force'][:,members]/z['Fmax'][members],0.).mean(axis=1)
    meta.update(changed_muscles=MUSCLES, native_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                parameters_sha256=hashlib.sha256(args.osim.read_bytes()).hexdigest(),
                maximum_equilibrium_error_N=max_error, maximum_fiber_step_s=args.max_step,
                residual_mean=float(z['residual'].mean()), residual_p95=float(np.percentile(z['residual'],95)),
                residual_scope='Bounded recruitment after changing four muscle force equations, same target torques.',
                maximum_optimality=float(max(optimality)),
                force_fields='Force is evaluated at the fitted activation. For the four nonlinear muscles, passive_force is the force at minimum activation and unit_force is the difference between full and minimum activation at the same incoming fiber state.')
    np.savez_compressed(out/'targets.npz',**z)
    (out/'metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf8')
    print(json.dumps({k:meta[k] for k in ('changed_muscles','residual_mean','residual_p95','maximum_equilibrium_error_N')}))


if __name__ == '__main__':
    main()

