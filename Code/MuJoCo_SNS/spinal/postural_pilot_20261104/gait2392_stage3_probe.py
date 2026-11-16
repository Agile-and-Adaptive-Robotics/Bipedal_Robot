"""Matched Stage 3 standing probe with pelvis and lumbar premotor corrections."""
import argparse
import json
import os
from pathlib import Path
import runpy
import sys

for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
HERE = Path(__file__).resolve().parent
SPINAL = HERE.parent
sys.path[:0] = [str(SPINAL / 'hunt_training'), str(SPINAL), str(HERE)]

import mujoco
import numpy as np
from scipy.optimize import lsq_linear

from model import build


AXES = ('pelvis_tilt', 'pelvis_list', 'pelvis_rotation',
        'lumbar_extension', 'lumbar_bending', 'lumbar_rotation')


class PremotorBalance:
    def __init__(self, network, muscles, model, data, parameters, dt, gain):
        self.modules = [build(backend='compact', parameters=parameters, dt=dt) for _ in AXES]
        self.ports = np.array([network.ports['POST_' + muscle] for muscle in muscles])
        self.actuators = np.array([model.actuator(muscle).id for muscle in muscles])
        self.qids = np.array([model.joint(name).qposadr[0] for name in AXES])
        self.vids = np.array([model.joint(name).dofadr[0] for name in AXES])
        self.goal = data.qpos[self.qids].copy()
        self.gain = gain
        probe = mujoco.MjData(model)
        probe.qpos[:] = data.qpos
        probe.act[:] = .08
        mujoco.mj_forward(model, probe)
        self.length = probe.actuator_length[self.actuators].copy()
        self.force_scale = np.maximum(model.actuator_gainprm[self.actuators, 2], 5.)
        baseline = probe.qacc[self.vids].copy()
        response = np.empty((len(AXES), len(muscles)))
        for j, actuator in enumerate(self.actuators):
            probe.act[actuator] += .01
            mujoco.mj_forward(model, probe)
            response[:, j] = (probe.qacc[self.vids] - baseline) / .01
            probe.act[actuator] -= .01
        self.basis = np.zeros((len(AXES), 2, len(muscles)))
        trunk = np.array([name.startswith(('ercspn_', 'intobl_', 'extobl_')) for name in muscles])
        for j in range(len(AXES)):
            eligible = trunk if j >= 3 else ~trunk
            local = response[:, eligible]
            design = np.vstack((local, 2. * np.eye(local.shape[1])))
            for k, sign in enumerate((1., -1.)):
                target = np.zeros(len(AXES))
                target[j] = 3. * sign
                solution = lsq_linear(design, np.r_[target, np.zeros(local.shape[1])],
                                      bounds=(0., .22), method='bvls')
                self.basis[j, k, eligible] = solution.x
        self.recruitment = self.basis / np.maximum(self.basis.sum(axis=2, keepdims=True), 1e-9)
        self.delay = np.tile(np.r_[self.goal, np.zeros(len(AXES))], (round(.032 / dt) + 1, 1))
        self.index = 0
        self.trace = []
        self.evidence = dict(axes=AXES, muscle_names=muscles, recruitment=self.basis.tolist(),
                             acceleration_response=np.einsum('im,jkm->jki', response, self.basis).tolist())

    def update(self, inputs, data):
        self.delay[self.index, :6] = data.qpos[self.qids]
        self.delay[self.index, 6:] = data.qvel[self.vids]
        self.index = (self.index + 1) % len(self.delay)
        pose, velocity = self.delay[self.index].reshape(2, 6)
        command = np.clip(-self.gain * (8. * (pose + .032 * velocity - self.goal)
                                        + 2. * velocity), -1., 1.)
        length = np.clip((data.actuator_length[self.actuators] - self.length) / .02, 0., 1.)
        speed = np.clip(data.actuator_velocity[self.actuators] / .1, 0., 1.)
        force = np.clip(-data.actuator_force[self.actuators] / self.force_scale, 0., 1.)
        motor = np.empty((6, 2))
        for j, module in enumerate(self.modules):
            module.u[:2] = module.command(np.maximum((command[j], -command[j]), 0.))
            module.u[2:4] = 1.6 * (self.recruitment[j] @ speed)
            module.u[4:6] = 1.2 * (self.recruitment[j] @ length)
            module.u[6:8] = 1.5 * (self.recruitment[j] @ force)
            module.compiled.forward(module.u)
            motor[j] = np.clip(module.compiled.V[[module.idx['MN_POS'], module.idx['MN_NEG']]] / 5., 0., 1.)
        correction = 5. * np.einsum('jk,jkm->m', motor, self.basis)
        inputs[self.ports] += np.minimum(correction, 1.5)
        if len(self.trace) < 4100:
            self.trace.append(np.r_[data.time, data.qpos[self.qids], command, motor.ravel(),
                                    np.linalg.norm(correction)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--fit', type=Path, required=True)
    parser.add_argument('--gain', type=float, default=0.)
    parser.add_argument('--support', type=float, default=.45)
    parser.add_argument('--spring-scale', type=float, default=1.)
    parser.add_argument('--push', type=float, default=30.)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    params = json.loads(args.fit.read_text())['parameters']
    params = {key: val for key, val in params.items() if key not in ('kp', 'kd')}
    import controller
    original = controller.build
    balance = None

    def factory(acts, dt=.002, **kwargs):
        network = original(acts, dt=dt, **kwargs)
        if args.gain == 0.:
            return network
        step = network.step

        def advance(u):
            nonlocal balance
            if network.live is not None:
                model, data = network.live
                if balance is None:
                    balance = PremotorBalance(network, list(acts), model, data, params, dt, args.gain)
                balance.update(u, data)
            return step(u)

        network.step = advance
        return network

    cfg = dict(stage=3, backend='compact', standing_solver='contact_symmetric', restriction=.25,
               support_fraction=args.support, vertical_stiffness_scale=args.spring_scale,
               start_hip_deg=-7., start_knee_deg=-12.,
               start_ankle_deg=0., start_lumbar_extension_deg=0., regional_contacts=True,
               push_N=args.push, lateral_push_N=args.push, push_time=4., timing_afferents=False,
               verification_capture=str(args.output / 'capture'), premotor_gain=args.gain,
               premotor_checkpoint=str(args.fit))
    config = args.output / 'body.json'
    config.write_text(json.dumps(cfg, indent=2))
    controller.build = factory
    step = mujoco.mj_step

    def guarded(model, data, *pos, **kw):
        for side in ('l', 'r'):
            joint = model.joint('knee_angle_' + side)
            angle = float(data.qpos[joint.qposadr[0]])
            velocity = float(data.qvel[joint.dofadr[0]])
            if velocity > 0. and angle + model.opt.timestep * velocity >= np.radians(-.1):
                (args.output / 'KNEE_GUARD.json').write_text(json.dumps(
                    dict(side=side, time_s=data.time, angle_deg=float(np.degrees(angle)))))
                raise RuntimeError('Predictive knee guard')
        return step(model, data, *pos, **kw)

    mujoco.mj_step = guarded
    sys.argv = [str(SPINAL / 'hunt_training' / 'physical_eval.py'), str(config)]
    try:
        runpy.run_path(sys.argv[0], run_name='__main__', init_globals={'print': lambda *a, **k: None})
    finally:
        controller.build = original
        mujoco.mj_step = step
        if balance is not None:
            np.savez_compressed(args.output / 'PREMOTOR_TRACE.npz', samples=balance.trace)
            (args.output / 'RECRUITMENT.json').write_text(json.dumps(balance.evidence, indent=2))


if __name__ == '__main__':
    main()
