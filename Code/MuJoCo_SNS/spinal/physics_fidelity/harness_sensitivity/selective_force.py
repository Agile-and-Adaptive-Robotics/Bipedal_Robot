"""Apply the selective muscle substitution within one frozen evaluator process."""
import hashlib
import json
from pathlib import Path
import runpy
import sys
import importlib.util

import mujoco
import numpy as np


FIDELITY = Path(__file__).resolve().parents[1]
REPO = FIDELITY.parents[3]
MODEL = REPO / 'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim'
MUSCLES = ('rect_fem_r', 'rect_fem_l', 'bifemsh_r', 'bifemsh_l')
ARMS = {'native': MUSCLES, 'selective': MUSCLES,
        'rect_fem': MUSCLES[:2], 'bifemsh': MUSCLES[2:]}


def array_digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def record_start(model, data, ids, muscles, path):
    others = np.setdiff1d(np.arange(model.nu), ids)
    start = {field: array_digest(getattr(data, field))
             for field in ('qpos', 'qvel', 'act', 'ctrl', 'xfrc_applied')}
    start.update({field: array_digest(getattr(model, field)[others])
                  for field in ('actuator_gainprm', 'actuator_biasprm',
                                'actuator_gaintype', 'actuator_biastype')})
    start.update(model_time=float(data.time), muscle_ids=ids.tolist(),
                 muscles=muscles, physics_dt=float(model.opt.timestep))
    path.write_text(json.dumps(start, indent=2), encoding='utf-8')


class SelectiveMuscleStep:
    """Advance the selected bilateral fiber states once per accepted MuJoCo step."""

    def __init__(self, model, data, ids, muscles, real_step):
        from thelen_dynamic import DynamicThelenBank, parameters_from_osim

        self.ids = ids
        self.real_step = real_step
        self.dt = float(model.opt.timestep)
        self.parameters = parameters_from_osim(MODEL, muscles, 0.0)
        self.bank = DynamicThelenBank(self.parameters, self.dt)
        self.geometry = mujoco.MjData(model)
        self.update_geometry(model, data)
        self.bank.initialize(self.geometry.actuator_length[ids].copy(), data.act[ids].copy())
        if not np.all(np.isfinite(self.bank.force)):
            raise RuntimeError('Nonfinite initial muscle force')
        model.actuator_gaintype[ids] = mujoco.mjtGain.mjGAIN_FIXED
        model.actuator_gainprm[ids, 0] = 0.0
        model.actuator_biastype[ids] = mujoco.mjtBias.mjBIAS_AFFINE
        model.actuator_biasprm[ids, 0] = -self.bank.force
        model.actuator_biasprm[ids, 1:3] = 0.0
        self.samples = []
        self.steps = 0

    def update_geometry(self, model, data):
        self.geometry.qpos[:] = data.qpos
        self.geometry.qvel[:] = data.qvel
        mujoco.mj_fwdPosition(model, self.geometry)
        mujoco.mj_fwdVelocity(model, self.geometry)

    def advance(self, model, data, *args, **kwargs):
        from thelen_dynamic import fiber_components

        ids = self.ids
        self.update_geometry(model, data)
        length = self.geometry.actuator_length[ids].copy()
        speed = self.geometry.actuator_velocity[ids].copy()
        activation = data.act[ids].copy()
        initial_state = self.bank.state.copy()
        force = self.bank.evaluate(length + self.dt * speed, activation).copy()
        slope = -self.dt * self.bank.tangent.copy()
        # Tensile force runs in the negative transmission direction.
        model.actuator_biasprm[ids, 0] = -force - slope * speed
        model.actuator_biasprm[ids, 2] = slope
        self.real_step(model, data, *args, **kwargs)
        self.update_geometry(model, data)
        self.bank.commit(self.geometry.actuator_length[ids].copy(), activation)
        if (not np.all(np.isfinite(data.qpos)) or not np.all(np.isfinite(self.bank.force))
                or np.max(self.bank.force / self.parameters[:, 0]) > 20):
            raise RuntimeError(f'Invalid muscle or body state at {data.time:g} s')
        if np.any(self.bank.bounds) or np.max(np.abs(self.bank.error)) > 1e-5:
            raise RuntimeError(f'Fiber equilibrium or length limit at {data.time:g} s')
        self.steps += 1
        if self.steps % 5 == 0:
            parts = np.empty((3, len(ids)))
            fiber_components(self.bank.state, initial_state, activation,
                             self.parameters, self.dt, parts)
            # Applied and committed force use opposite signs in MuJoCo's convention.
            self.samples.append(np.r_[data.time, data.actuator_force[ids], self.bank.force,
                                      parts.ravel(), self.bank.error, self.bank.bounds])


def main():
    arm, script, config = sys.argv[1:4]
    corrected = sys.argv[4:] == ['--corrected-start']
    if not corrected and sys.argv[4:]:
        raise ValueError(sys.argv[4:])
    if arm not in ARMS:
        raise ValueError(arm)
    muscles = ARMS[arm]
    evaluator = Path(script)
    spinal = evaluator.parent.parent
    sys.path[:0] = [str(evaluator.parent), str(spinal), str(spinal / 'functional_gait2392')]
    # Load this experiment's checked implementation before the frozen evaluator
    # adds its own historical physics_fidelity modules to the search path.
    sys.path.insert(0, str(FIDELITY))
    import thelen_dynamic
    if Path(thelen_dynamic.__file__).resolve() != FIDELITY / 'thelen_dynamic.py':
        raise RuntimeError('Incorrect Thelen force-law source')
    if corrected:
        source = FIDELITY.parent / 'joint_projection.py'
        spec = importlib.util.spec_from_file_location('source_joint_projection', source)
        projection = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(projection)
        import runner
        runner._project_followers = projection.project_joint_followers
    target = Path(config).parent
    real_step = mujoco.mj_step
    started = False
    fiber_step = None

    def intercept(model, data, *args, **kwargs):
        nonlocal started, fiber_step
        # The archived evaluator also calls MuJoCo when preparing the apparatus.
        # Change forces only in runner.main's coupled closed-loop integration.
        frame = sys._getframe(1)
        if frame.f_code.co_name != 'main' or Path(frame.f_code.co_filename).name != 'runner.py':
            return real_step(model, data, *args, **kwargs)
        if not started:
            ids = np.array([model.actuator(name).id for name in muscles])
            verified = (model.nu == model.na == 92
                        and np.all(model.actuator_gaintype[ids] == mujoco.mjtGain.mjGAIN_MUSCLE)
                        and np.all(model.actuator_biastype[ids] == mujoco.mjtBias.mjBIAS_MUSCLE)
                        and np.all(model.actuator_dynprm[ids, 0] > 0))
            if not verified or (hasattr(model, 'actuator_actadr')
                                and not np.array_equal(model.actuator_actadr[ids], ids)):
                raise ValueError('Expected the archived 92-muscle native Gait2392 mechanism')
            record_start(model, data, ids, muscles, target / 'controlled_start.json')
            if arm != 'native':
                fiber_step = SelectiveMuscleStep(model, data, ids, muscles, real_step)
            started = True
        if fiber_step is not None:
            return fiber_step.advance(model, data, *args, **kwargs)
        return real_step(model, data, *args, **kwargs)

    mujoco.mj_step = intercept
    sys.argv = [script, config]
    try:
        runpy.run_path(script, run_name='__main__')
    finally:
        mujoco.mj_step = real_step
        if fiber_step is not None:
            np.savez_compressed(target / 'selective_force_trace.npz',
                                t_applied_corrected_active_passive_damping_equilibrium_bounds=
                                np.array(fiber_step.samples), names=muscles)


if __name__ == '__main__':
    main()
