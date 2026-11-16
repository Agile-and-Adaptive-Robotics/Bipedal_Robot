"""Fixed-base sagittal serial balance with a torque reference and SNS transfer."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)

import mujoco
import numpy as np
from scipy.linalg import expm, solve_continuous_are

from model import build


DT = .001
CAPACITY = 160.
FILTER = .03


def model_xml(count, axes=1):
    links = []
    actuators = []
    for j, (length, mass) in enumerate(zip((.45, .42, .35), (4., 5., 10.))):
        if j == count:
            break
        offset = (0., .45, .42)[j]
        joints = []
        for axis in ('z', 'x')[:axes]:
            joint = f'joint{j}_{axis}'
            direction = '0 0 1' if axis == 'z' else '1 0 0'
            joints.append(f'<joint name="{joint}" axis="{direction}" damping=".2"/>')
            for sign in (1, -1):
                actuators.append(f'<general joint="{joint}" gear="{sign}" dyntype="filter" '
                                 f'dynprm="{FILTER}" gainprm="{CAPACITY}" ctrlrange="0 1" ctrllimited="true"/>')
        links.append(f'<body name="link{j}" pos="0 {offset} 0">{"".join(joints)}'
                     f'<geom type="capsule" fromto="0 0 0 0 {length} 0" '
                     f'size=".035" mass="{mass}" contype="0" conaffinity="0"/>')
    return (f'<mujoco model="serial{count}"><option timestep="{DT}" gravity="0 -9.81 0" '
            f'integrator="implicitfast"/><worldbody>{"".join(links)}'
            f'{"</body>" * count}</worldbody><actuator>{"".join(actuators)}</actuator></mujoco>')


def linear_model(model):
    n = model.nv
    data = mujoco.MjData(model)
    a = np.zeros((3 * n, 3 * n))
    a[:n, n:2 * n] = np.eye(n)
    a[2 * n:, 2 * n:] = -np.eye(n) / FILTER
    eps = 1e-5
    for i in range(3 * n):
        target = (data.qpos, data.qvel, data.qfrc_applied)[i // n]
        target[i % n] = eps
        mujoco.mj_forward(model, data)
        plus = data.qacc.copy()
        target[i % n] = -eps
        mujoco.mj_forward(model, data)
        a[n:2 * n, i] = (plus - data.qacc) / (2 * eps)
        target[i % n] = 0.
    b = np.zeros((3 * n, n))
    b[2 * n:] = np.eye(n) / FILTER
    return a, b


def feedback_gain(model, effort, neural_tau=None, distal_effort=1., velocity_weight=1.):
    a, b = linear_model(model)
    n = model.nv
    weights = [180., 10. * velocity_weight, .002]
    if neural_tau is not None:
        extended = np.zeros((5 * n, 5 * n))
        extended[:3 * n, :3 * n] = a
        extended[:3 * n, 3 * n:4 * n] = b
        extended[3 * n:4 * n, 3 * n:4 * n] = -np.eye(n) / .01
        extended[3 * n:4 * n, 4 * n:] = np.eye(n) / .01
        extended[4 * n:, 4 * n:] = -np.eye(n) / neural_tau
        a = extended
        b = np.zeros((5 * n, n))
        b[4 * n:] = np.eye(n) / neural_tau
        weights += [.002, .002]
    q = np.diag(np.repeat(weights, n))
    r = np.eye(n) * effort
    r[-1, -1] *= distal_effort
    gain = np.linalg.solve(r, b.T @ solve_continuous_are(a, b, q, r))
    return gain


def evaluate(model, gain, mode, sign, push, neural, delay_s=0., prediction=0., initial_scale=1.,
             push_axis='last', model_predictor=False, disturbance_gain=0., capture=None):
    n = model.nv
    data = mujoco.MjData(model)
    data.qpos[:] = sign * initial_scale * np.linspace(.028, -.015, n)
    modules = [build(backend='compact', parameters=neural) for _ in range(n)] if mode.startswith('SNS') else []
    mn = [[module.idx['MN_POS'], module.idx['MN_NEG']] for module in modules]
    pp = [[module.idx['PP_POS'], module.idx['PP_NEG']] for module in modules]
    activation = data.act.reshape(n, 2)
    control = data.ctrl.reshape(n, 2)
    state = np.zeros(gain.shape[1])
    delay_steps = round(delay_s / DT)
    sensors = np.zeros((delay_steps + 1, 2 * n))
    sensors[:, :n] = data.qpos
    use_predictor = bool(delay_steps and model_predictor)
    if use_predictor:
        a, _ = linear_model(model)
        augmented = np.zeros((3 * n, 3 * n))
        augmented[:2 * n] = a[:2 * n]
        discretized = expm(augmented * DT)
        transition = discretized[:2 * n, :2 * n]
        torque_response = discretized[:2 * n, 2 * n:]
        weights = [np.linalg.matrix_power(transition, delay_steps - 1 - i) @ torque_response
                   for i in range(delay_steps)]
        lift = np.column_stack(weights)
        horizon = np.linalg.matrix_power(transition, delay_steps)
        history = np.zeros((2 * delay_steps, n))
        predicted = np.empty(2 * n)
        correction = np.empty(2 * n)
        if disturbance_gain:
            # Integrate a constant estimated disturbance acceleration over the delay.
            augmented[:2 * n, 2 * n:] = np.vstack((np.zeros((n, n)), np.eye(n)))
            response = expm(augmented * DT)[:2 * n, 2 * n:]
            disturbance_horizon = sum(np.linalg.matrix_power(transition, i) @ response
                                      for i in range(delay_steps))
            residual = np.zeros(n)
            previous_velocity = np.zeros(n)
            residual_alpha = -np.expm1(-DT / .008)
    direction = np.array([1., -1.])
    push_vector = np.zeros(n)
    if push_axis == 'diagonal':
        push_vector[-2:] = sign * push / np.sqrt(2.)
    else:
        push_vector[-2 if push_axis == 'other' else -1] = sign * push
    log = np.empty((1200, 1 + 2 * n))
    hold = 0
    settled = None
    peak = 0.
    used = 0
    started = time.perf_counter()
    for step in range(6000):
        sensors[step % len(sensors), :n] = data.qpos
        sensors[step % len(sensors), n:] = data.qvel
        state[:2 * n] = sensors[(step - delay_steps) % len(sensors)]
        if use_predictor:
            index = step % delay_steps
            if step >= delay_steps:
                np.matmul(horizon, state[:2 * n], out=predicted)
                np.matmul(lift, history[index:index + delay_steps].reshape(-1), out=correction)
                if disturbance_gain:
                    velocity = state[n:2 * n]
                    expected = a[n:2 * n, :2 * n] @ state[:2 * n] + a[n:2 * n, 2 * n:3 * n] @ history[index]
                    measured = (velocity - previous_velocity) / DT
                    residual += residual_alpha * (measured - expected - residual)
                    previous_velocity[:] = velocity
                    correction += disturbance_gain * (disturbance_horizon @ residual)
                state[:2 * n] += prediction * (predicted + correction - state[:2 * n])
        elif delay_steps and prediction:
            state[:n] += prediction * delay_s * state[n:2 * n]
        state[2 * n:3 * n] = CAPACITY * (activation[:, 0] - activation[:, 1])
        if use_predictor:
            history[index] = history[index + delay_steps] = state[2 * n:3 * n]
        if len(state) == 5 * n:
            for j, module in enumerate(modules):
                motor = np.clip(module.compiled.V[mn[j]] / 5., 0., 1.)
                premotor = np.clip(module.compiled.V[pp[j]] / 5., 0., 1.) * (5. / 3.)
                premotor = 1.6 * premotor / (1. + premotor)
                state[3 * n + j] = .95 * CAPACITY * (motor[0] - motor[1])
                state[4 * n + j] = .95 * CAPACITY * (premotor[0] - premotor[1])
        torque = -gain @ state
        if mode == 'diagonal PD':
            torque = -np.diag(gain[:, :n]) * data.qpos - np.diag(gain[:, n:2 * n]) * data.qvel
        desired = np.clip(np.column_stack((torque, -torque)) / CAPACITY, 0., 1.)
        if modules:
            for j, module in enumerate(modules):
                module.u[:2] = module.command(desired[j])
                # Signed virtual actuator length from constant-radius joint transmission.
                module.u[2:4] = 1.6 * np.clip(-direction * data.qvel[j] * .05 / .1, 0., 1.)
                module.u[4:6] = 1.2 * np.clip(-direction * data.qpos[j] * .05 / .02, 0., 1.)
                module.u[6:8] = 1.5 * activation[j]
                module.compiled.forward(module.u)
                control[j] = .05 + .95 * np.clip(module.compiled.V[mn[j]] / 5., 0., 1.)
        else:
            control[:] = desired
        if step == 1500:
            data.qfrc_applied[:] = push_vector
        elif step == 1550:
            data.qfrc_applied[:] = 0.
        mujoco.mj_step(model, data)
        peak = max(peak, float(np.max(abs(data.qpos))))
        if not np.isfinite(data.qpos).all() or peak > .5:
            break
        if step % 5 == 0:
            log[used] = np.r_[data.time, data.qpos, data.qvel]
            if capture is not None:
                capture['actuator_command'][used] = control
                capture['actuator_activation'][used] = activation
                capture['push_Nm'][used] = data.qfrc_applied
                for j, module in enumerate(modules):
                    capture['premotor_mV'][used, j] = module.compiled.V[pp[j]]
                    capture['motor_mV'][used, j] = module.compiled.V[mn[j]]
            used += 1
        if step >= 1550:
            in_band = (np.max(abs(data.qpos)) < np.radians(1.)
                       and np.max(abs(data.qvel)) < np.radians(5.))
            hold = hold + 1 if in_band else 0
            if hold == 500 and settled is None:
                settled = data.time - 1.55 - .499
    log = log[:used]
    tail = log[log[:, 0] >= 4.]
    rms = float(np.degrees(np.sqrt(np.mean(tail[:, 1:n + 1] ** 2)))) if len(tail) else None
    complete = data.time >= 6. - DT / 2
    return dict(links=model.nbody-1, dof=n, mode=mode, mirrored_sign=sign,
                push_Nm=push, push_axis=push_axis, delay_s=delay_s, complete=complete,
                peak_angle_deg=float(np.degrees(peak)), final_rms_deg=rms,
                settling_s=settled, passed=bool(complete and rms is not None and rms < 1.
                                              and settled is not None and settled < 3.),
                simulated_s=float(data.time), wall_s=time.perf_counter() - started), log


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--single-fit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--lag-sweep', action='store_true')
    parser.add_argument('--axes', type=int, choices=(1, 2), default=1)
    parser.add_argument('--planar-fit-dir', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    fit = json.loads(args.single_fit.read_text())['parameters']
    neural = {key: val for key, val in fit.items() if key not in ('kp', 'kd')}
    results = []
    for count in (2, 3):
        xml = model_xml(count, args.axes)
        (args.output / f'links{count}.xml').write_text(xml)
        model = mujoco.MjModel.from_xml_string(xml)
        gain = feedback_gain(model, .025)
        if args.axes == 2 and args.planar_fit_dir:
            prior = json.loads((args.planar_fit_dir / f'links{count}_selected.json').read_text())
            lag_gain = feedback_gain(model, prior['effort'], prior['modeled_neural_tau_s'])
            for mode, controller in (('full-state LQR', gain), ('SNS lag-aware', lag_gain)):
                for axis in ('last', 'other', 'diagonal'):
                    for sign in (1, -1):
                        row, trace = evaluate(model, controller, mode, sign, 8., neural, push_axis=axis)
                        results.append(row)
                        np.savez_compressed(args.output / f'links{count}_{mode}_{axis}_{sign}.npz',
                                            time_state=trace)
            (args.output / 'RESULTS.json').write_text(json.dumps(results, indent=2))
            (args.output / f'links{count}_gain.json').write_text(json.dumps(lag_gain.tolist()))
            continue
        if args.lag_sweep:
            best = None
            for effort in (.01, .025, .07, .2, .6, 2.):
                for tau in (.008, .015, .025):
                    candidate = feedback_gain(model, effort, neural_tau=tau)
                    row, _ = evaluate(model, candidate, 'SNS lag-aware', 1, 8., neural)
                    row.update(effort=effort, modeled_neural_tau_s=tau)
                    results.append(row)
                    score = (not row['passed'], -(row['simulated_s']),
                             row['settling_s'] if row['settling_s'] is not None else 1e3,
                             row['peak_angle_deg'])
                    if best is None or score < best[0]:
                        best = (score, row, candidate)
            _, selected, gain = best
            (args.output / f'links{count}_selected.json').write_text(json.dumps(selected, indent=2))
            for sign in (1, -1):
                for push, delay in ((8., 0.), (12., .01), (16., .02)):
                    row, trace = evaluate(model, gain, 'SNS lag-aware', sign, push, neural, delay)
                    results.append(row | {'validation': True})
                    np.savez_compressed(args.output / f'links{count}_validation_{sign}_{push:g}.npz',
                                        time_state=trace)
            (args.output / 'RESULTS.json').write_text(json.dumps(results, indent=2))
            (args.output / f'links{count}_gain.json').write_text(json.dumps(gain.tolist()))
            continue
        for mode in ('diagonal PD', 'full-state LQR', 'SNS'):
            for sign in (1, -1):
                row, trace = evaluate(model, gain, mode, sign, 8., neural)
                results.append(row)
                np.savez_compressed(args.output / f'links{count}_{mode.replace(" ", "_")}_{sign}.npz',
                                    time_state=trace)
        (args.output / 'RESULTS.json').write_text(json.dumps(results, indent=2))
        (args.output / f'links{count}_gain.json').write_text(json.dumps(gain.tolist()))
    print(json.dumps(results))


if __name__ == '__main__':
    main()
