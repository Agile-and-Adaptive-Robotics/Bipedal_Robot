"""Model-based parity, loaded FK, and independent velocity checks for desired IK."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import mujoco
from scipy.spatial.transform import Rotation
from supported_step_teacher import NativeCartesianTeacher


def validate(model_folder, start, end, spec_path, legacy_path, replay, out):
    spec = json.loads(spec_path.read_text())
    frames = []
    for folder in (start, end, start):
        with np.load(folder/'source_recruitment.npz') as z:
            frames.append({'qpos': z['qpos'].copy()})
    module_spec = importlib.util.spec_from_file_location('legacy_cartesian', legacy_path)
    legacy = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(legacy)
    m = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
    d = mujoco.MjData(m)
    with np.load(replay/'transition_sparse.npz') as z:
        samples = z['samples'].copy()
    times = samples[:, 0]
    positions = samples[:, 1:1+m.nq]
    velocities = np.gradient(positions, times, axis=0)
    parity_q = parity_v = 0.
    for hybrid in (False, True):
        old = legacy.NativeCartesianTeacher(m, frames, spec, hybrid)
        new = NativeCartesianTeacher(m, frames, spec, hybrid)
        for k in np.linspace(0, len(samples)-1, 21).astype(int):
            t = np.clip(times[k]-.5, 0., 4.)
            index = min(int(t/2.), 1)
            phase = np.clip((t-2.*index)/2., 0., 1.)
            blend = phase**3*(10.+phase*(-15.+6.*phase))
            rate = 15.*phase**2*(1.-phase)**2
            q = frames[index]['qpos']+blend*(frames[index+1]['qpos']-frames[index]['qpos'])
            v = rate*(frames[index+1]['qpos']-frames[index]['qpos'])
            d.qpos[:] = positions[k]; d.qvel[:] = velocities[k]; d.time = times[k]
            mujoco.mj_fwdPosition(m, d)
            qo, vo = old.solve(d, index, blend, rate, q, v)
            qn, vn = new.solve(d, index, blend, rate, q, v)
            parity_q = max(parity_q, float(np.max(abs(qn-qo))))
            parity_v = max(parity_v, float(np.max(abs(vn-vo))))
    teacher = NativeCartesianTeacher(m, frames, spec, preserve_loaded_followers=True)
    keyframe_error = 0.
    for k, frame in enumerate(frames):
        d.qpos[:] = frame['qpos']; mujoco.mj_fwdPosition(m, d)
        keyframe_error = max(keyframe_error, float(np.max(abs(d.xpos[teacher.feet]-teacher.position[k]))))
    d.qpos[:] = frames[0]['qpos']; mujoco.mj_fwdPosition(m, d)
    teacher.arms.update(d)
    d.qvel[:] = 0.; d.qvel[:6] = np.arange(1, 7)*1e-4
    d.qvel[teacher.arms.rows] = np.sin(np.arange(len(teacher.arms.rows)))*1e-3
    deformation_rate = np.sin(np.arange(len(teacher.arms.dependent)))*1e-5
    d.qvel[teacher.arms.dependent] = teacher.arms.E[teacher.arms.dependent]@d.qvel[teacher.arms.rows]+deformation_rate
    saved_q = d.qpos.copy(); saved_v = d.qvel.copy()
    q, v = teacher.solve(d, 0, 0., 0., frames[0]['qpos'], np.zeros(m.nv))
    unchanged = bool(np.array_equal(d.qpos, saved_q) and np.array_equal(d.qvel, saved_v))
    v[:6] = d.qvel[:6]
    probe = mujoco.MjData(m); probe.qpos[:] = q; mujoco.mj_fwdPosition(m, probe)
    p0 = probe.xpos[teacher.feet].copy(); r0 = probe.xmat[teacher.feet].reshape(2, 3, 3).copy()
    epsilon = 1e-4
    probe.qpos[:] = q+epsilon*v; mujoco.mj_fwdPosition(m, probe)
    linear_error = float(np.max(abs((probe.xpos[teacher.feet]-p0)/epsilon)))
    angular_error = max(float(np.linalg.norm(Rotation.from_matrix(probe.xmat[b].reshape(3, 3)@r0[k].T).as_rotvec()/epsilon)) for k, b in enumerate(teacher.feet))
    profile = NativeCartesianTeacher(m, frames, spec, preserve_loaded_followers=True, tangential_clearance=.005)
    profile_endpoint = profile_derivative = profile_loaded_speed = 0.
    for index in range(len(frames)-1):
        for side in range(2):
            for blend in np.linspace(0., 1., 101):
                position, orientation, velocity, gap = profile.target_kinematics(index, side, blend, 1.)
                if blend in (0., 1.):
                    k = index+int(blend)
                    profile_endpoint = max(profile_endpoint, float(np.max(abs(position-profile.position[k][side]))), float(np.max(abs(orientation-profile.rotation[k][side]))))
                else:
                    h = 1e-6
                    plus = profile.target_kinematics(index, side, blend+h, 1.)
                    minus = profile.target_kinematics(index, side, blend-h, 1.)
                    measured = np.r_[(plus[0]-minus[0])/(2*h), Rotation.from_matrix(plus[1]@minus[1].T).as_rotvec()/(2*h)]
                    profile_derivative = max(profile_derivative, float(np.max(abs(measured-velocity))))
                if gap <= .001:
                    profile_loaded_speed = max(profile_loaded_speed, float(np.linalg.norm(velocity[:2])))
    result = dict(default_off_replay_qpos_max_error=parity_q, default_off_replay_qvel_max_error=parity_v,
                  loaded_keyframe_world_position_error_m=keyframe_error,
                  independent_stationary_foot_velocity_error_mps=linear_error,
                  independent_stationary_foot_angular_velocity_error_radps=angular_error,
                  live_state_unchanged=unchanged, replay_samples=42,
                  profile_endpoint_error=profile_endpoint, profile_independent_derivative_error=profile_derivative,
                  profile_low_clearance_tangent_per_unit_blend=profile_loaded_speed,
                  legacy_code=str(legacy_path), replay=str(replay), loaded_diagnostics=teacher.diagnostics())
    result['accepted'] = bool(parity_q == 0. and parity_v == 0. and keyframe_error == 0. and unchanged and linear_error < 1e-7 and angular_error < 1e-7 and profile_endpoint < 1e-12 and profile_derivative < 1e-7 and profile_loaded_speed < .005)
    out.mkdir(parents=True, exist_ok=True)
    target = out/'LOADED_CARTESIAN_VALIDATION.json'
    if target.exists():
        raise ValueError('Preserve existing validation')
    target.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if not result['accepted']:
        raise SystemExit(1)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('model', 'start', 'end', 'spec', 'legacy', 'replay', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    validate(a.model, a.start, a.end, a.spec, a.legacy, a.replay, a.out)
