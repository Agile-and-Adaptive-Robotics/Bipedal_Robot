"""Offline native-activation/source-fiber prediction on immutable trajectories.

The oracle replay validates implementation. Constant-command, linear-length
forecasts are causal approximations, not accepted inverse-motion references.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent/'functional_gait2392')]
import numpy as np
import mujoco
from numba import njit
from mechanics import CoupledArms, JOINTS
from thelen_dynamic import DynamicThelenBank


@njit(cache=True)
def activation_step(activation, excitation, parameters, dt):
    for i in range(len(activation)):
        a = activation[i]
        scale = .5+1.5*a
        tau = parameters[i, 0]*scale if excitation[i] > a else parameters[i, 1]/scale
        activation[i] += dt*(excitation[i]-a)/tau


class SourcePredictor:
    """Preallocated forecast with original native activation and source fibers."""
    def __init__(self, parameters, activation_parameters, dt, filter_tau):
        if np.any(activation_parameters[:, 2] != 0.):
            raise ValueError('Only the captured zero-smoothing native activation law is validated')
        self.dt = dt
        self.parameters = parameters
        self.activation_parameters = np.ascontiguousarray(activation_parameters[:, :3])
        self.beta = 1.-np.exp(-dt/filter_tau) if filter_tau else 1.
        self.bank = DynamicThelenBank(parameters, dt)
        self.activation = np.zeros(len(parameters))
        self.excitation = np.zeros_like(self.activation)
        self.length = np.zeros_like(self.activation)
        self.calls = 0
        self.invalid = 0
        self.negative_gains = 0
        self.max_residual = 0.

    def predict(self, state, activation, excitation, command, initial_length, velocity, steps, lengths=None, excitations=None):
        self.bank.state[:] = state
        self.activation[:] = activation
        self.excitation[:] = excitation
        for k in range(steps):
            if excitations is None:
                self.excitation += self.beta*(command-self.excitation)
            else:
                self.excitation[:] = excitations[k]
            if lengths is None:
                self.length[:] = initial_length+(k+1)*self.dt*velocity
            else:
                self.length[:] = lengths[k]
            self.bank.commit(self.length, self.activation)
            self.max_residual = max(self.max_residual, float(np.max(abs(self.bank.error))))
            if np.any(self.bank.bounds) or not np.all(np.isfinite(self.bank.force)) or np.max(abs(self.bank.error))>1e-5 or np.any(self.bank.force>20.*self.parameters[:,0]):
                return False
            activation_step(self.activation, self.excitation, self.activation_parameters, self.dt)
        return True

    def linearize_command(self, state, activation, excitation, command, length, velocity, steps, epsilon=.005):
        """Diagonal command derivative conditional on the causal length forecast."""
        self.calls += 1
        base = np.asarray(command).copy()
        lo = np.maximum(.01, base-epsilon)
        hi = np.minimum(.98, base+epsilon)
        if np.any(base < .01-1e-12) or np.any(base > .98+1e-12) or np.any(hi <= lo):
            self.invalid += 1
            return None
        valid = self.predict(state, activation, excitation, base, length, velocity, steps)
        f0 = self.bank.force.copy()
        valid = self.predict(state, activation, excitation, lo, length, velocity, steps) and valid
        fl = self.bank.force.copy()
        valid = self.predict(state, activation, excitation, hi, length, velocity, steps) and valid
        derivative = (self.bank.force-fl)/(hi-lo)
        tolerance = 1e-6*np.maximum(1., self.parameters[:, 0])
        self.negative_gains += int(np.any(derivative < -tolerance))
        valid = valid and np.all(np.isfinite(derivative)) and np.all(derivative >= -tolerance)
        if not valid:
            self.invalid += 1
            return None
        return f0, np.maximum(0., derivative), base

    def diagnostics(self):
        return dict(calls=self.calls, invalid_forecasts=self.invalid, negative_gain_forecasts=self.negative_gains,
                    maximum_constitutive_residual_N=self.max_residual,
                    scope='Independent conditional per-muscle forecasts; native activation,filter,and source fiber state. Linear length forecast,not future contact-state knowledge.')


def audit(model_folder, trace_folder, out, samples=101):
    if (out/'ACTUATOR_PREDICTION_AUDIT.json').exists():
        raise ValueError('Preserve existing prediction audit')
    code_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with np.load(trace_folder/'training_reference.npz') as archive:
        z = {k: archive[k] for k in archive.files}
    m = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
    names = z['muscles'].tolist()
    if names != [m.actuator(i).name for i in range(m.nu)]:
        raise ValueError('Actuator order mismatch')
    if not z['force_channel_enabled'].all() or np.any(m.actuator_dyntype != mujoco.mjtDyn.mjDYN_MUSCLE):
        raise ValueError('Explicit captured source92 muscle policy required')
    dt = float(np.median(np.diff(z['t'])))
    parameters = z['parameters']
    activation_parameters = m.actuator_dynprm[:, :3].copy()
    manifest = json.loads((trace_folder/'MANIFEST.json').read_text())
    filter_tau = float(manifest['source_feedback_allocator']['excitation_filter_tau_s'])
    predictor = SourcePredictor(parameters, activation_parameters, dt, filter_tau)
    grid_error = 0.
    for a in np.linspace(0., 1., 23):
        for u in np.linspace(0., 1., 23):
            predicted = np.array([a])
            activation_step(predicted, np.array([u]), activation_parameters[:1], dt)
            expected = a+dt*mujoco.mju_muscleDynamics(u, a, activation_parameters[0])
            grid_error = max(grid_error, float(abs(predicted[0]-expected)))
    if grid_error > 1e-12:
        raise ValueError('Compiled native activation formula failed primary API check')
    activation_error = 0.
    alpha = np.empty(m.nu)
    for k in range(len(z['t'])):
        alpha[:] = z['activation_for_force'][k]
        activation_step(alpha, z['excitation'][k], activation_parameters, dt)
        activation_error = max(activation_error, float(np.max(abs(alpha-z['activation'][k]))))
    alpha_alignment = float(np.max(abs(z['activation_for_force'][1:]-z['activation'][:-1])))
    if activation_error > 1e-10 or alpha_alignment > 1e-12:
        raise ValueError('Recorded/native activation timing mismatch')
    state = np.sqrt(np.maximum(0., z['fiber_length_m']**2-parameters[:, 3]**2))
    horizons = np.array([.01, .02, .05])
    steps = np.round(horizons/dt).astype(int)
    indices = np.linspace(max(1, round(.5/dt)), len(z['t'])-steps.max()-1, samples).astype(int)
    modes = ('actual_current', 'static_current', 'static_gain_command', 'oracle', 'held_command_oracle_length', 'held_command_linear_length')
    torque_error = np.full((len(modes), len(steps), len(indices), len(JOINTS)), np.nan)
    force_error = np.full((len(modes), len(steps), len(indices), m.nu), np.nan)
    invalid = np.zeros((len(modes), len(steps)), int)
    negative = np.zeros_like(invalid)
    command_span = np.empty((len(steps), len(indices)))
    contact_transition = np.zeros_like(command_span, dtype=bool)
    beta = predictor.beta
    target_commands = z['excitation'][1:]+(1./beta-1.)*(z['excitation'][1:]-z['excitation'][:-1])
    if np.max(target_commands) > .9800001 or np.min(target_commands) < .0099999:
        raise ValueError('Inferred causal engineering command exceeds original bounds')
    static = DynamicThelenBank(parameters, dt)
    d = mujoco.MjData(m)
    arms = CoupledArms(m)
    for j, k in enumerate(indices):
        length = z['mtu_length_m'][k]
        velocity = (length-z['mtu_length_m'][k-1])/dt
        command = target_commands[k-1]
        static.initialize(length, z['activation'][k])
        static_force = static.force.copy()
        static_valid = not np.any(static.bounds)
        lo = np.maximum(.01, z['activation'][k]-.005)
        hi = np.minimum(.98, z['activation'][k]+.005)
        static.initialize(length, lo)
        fl = static.force.copy()
        gain_valid = not np.any(static.bounds)
        static.initialize(length, hi)
        gain = np.maximum(0., (static.force-fl)/(hi-lo))
        gain_valid = gain_valid and not np.any(static.bounds)
        instantaneous = z['source_force_N'][k]+gain*(command-z['activation'][k])
        for hi, count in enumerate(steps):
            future = slice(k+1, k+count+1)
            actual = z['source_force_N'][k+count]
            d.qpos[:] = z['qpos'][k+count]
            mujoco.mj_fwdPosition(m, d)
            rho = arms.update(d).copy()
            command_span[hi, j] = np.max(abs(target_commands[k:k+count]-command))
            contact_transition[hi, j] = np.any((z['normal_side_BW'][k] <= .025) & (z['normal_side_BW'][k+count] > .025))
            for mode in range(len(modes)):
                if mode == 0:
                    force = z['source_force_N'][k]
                    valid = True
                elif mode == 1:
                    force = static_force
                    valid = static_valid
                elif mode == 2:
                    force = instantaneous
                    valid = gain_valid
                else:
                    valid = predictor.predict(state[k], z['activation'][k], z['excitation'][k], command,
                                              length, velocity, int(count),
                                              lengths=z['mtu_length_m'][future] if mode in (3, 4) else None,
                                              excitations=z['excitation'][future] if mode == 3 else None)
                    force = predictor.bank.force
                if not valid:
                    invalid[mode, hi] += 1
                    continue
                negative[mode, hi] += int(np.any(force < 0.))
                force_error[mode, hi, j] = force-actual
                torque_error[mode, hi, j] = rho@(force-actual)
    oracle_error = float(np.nanmax(abs(force_error[3])))
    if oracle_error > 1e-5:
        raise ValueError('Oracle replay failed exact source-force consistency')
    rows = []
    for mode, name in enumerate(modes):
        for hi, horizon in enumerate(horizons):
            stable = command_span[hi] <= .01
            rows.append(dict(mode=name, horizon_s=float(horizon), invalid_predictions=int(invalid[mode, hi]),negative_force_predictions=int(negative[mode, hi]),
                             normalized_force_RMS_F0=float(np.sqrt(np.nanmean((force_error[mode, hi]/parameters[:, 0])**2))),
                             joint_torque_RMS_Nm=np.sqrt(np.nanmean(torque_error[mode, hi]**2, axis=0)).tolist(),
                             stable_command_samples=int(stable.sum()),
                             stable_command_torque_RMS_Nm=float(np.sqrt(np.nanmean(torque_error[mode, hi, stable]**2))),
                             touchdown_samples=int(contact_transition[hi].sum()),
                             touchdown_torque_RMS_Nm=float(np.sqrt(np.nanmean(torque_error[mode, hi, contact_transition[hi]]**2))) if contact_transition[hi].any() else None))
    result = dict(source_trace=str(trace_folder), source_code_sha256=code_hash, samples=len(indices),
                  source_context='Historical mechanical failure,used only for constitutive/activation replay; not a current teacher reference.',
                  native_activation_grid_error=grid_error, native_activation_trace_error=activation_error,
                  native_activation_timing_error=alpha_alignment, oracle_force_max_error_N=oracle_error,
                  joint_names=JOINTS, rows=rows,
                  scope='Oracle uses recorded future lengths/excitation only to validate implementation. Held-command linear-length mode uses only current/past state. Future moment arms evaluate torque error,not the forecast. No physical model or controller is changed.')
    out.mkdir(parents=True, exist_ok=True)
    (out/'ACTUATOR_PREDICTION_AUDIT.json').write_text(json.dumps(result, indent=2))
    np.savez_compressed(out/'prediction_windows.npz',indices=indices,t=z['t'][indices],horizons=horizons,
                        modes=modes,force_error_N=force_error,torque_error_Nm=torque_error,
                        command_span=command_span,contact_transition=contact_transition)
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))


class AugmentedTransition:
    """Preallocated native rollout of q,v,activation,fibers,and filtering."""
    def __init__(self, model, parameters, filter_tau, max_steps=16,backend_class=None):
        import copy
        from thelen_dynamic import ConstitutiveBackend
        self.model = copy.copy(model)
        self.data = mujoco.MjData(self.model)
        mujoco.mj_forward(self.model, self.data)
        self.backend = (ConstitutiveBackend if backend_class is None else backend_class)(self.model, self.data, parameters)
        self.mask = mujoco.mjtState.mjSTATE_INTEGRATION
        self.saved_state = np.empty(mujoco.mj_stateSize(self.model, self.mask))
        self.saved_bias = self.model.actuator_biasprm.copy()
        self.saved_fiber = np.empty(self.model.nu)
        self.saved_force = np.empty(self.model.nu)
        self.saved_excitation = np.empty(self.model.nu)
        self.excitation = np.empty(self.model.nu)
        dt = float(self.model.opt.timestep)
        self.beta = 1.-np.exp(-dt/filter_tau) if filter_tau else 1.
        self.physics_steps = 0
        self.source_rom = None
        # Diagnostic buffers do not alter native activation or fiber branches.
        self.activation_distance = np.empty((max_steps, self.model.nu))
        self.fiber_velocity = np.empty((max_steps, self.model.nu))

    def load_snapshot(self, z, index):
        self.saved_state[:] = z['integration_state'][index]
        self.saved_bias[:] = z['actuator_biasprm'][index]
        self.saved_fiber[:] = z['fiber_state'][index]
        self.saved_force[:] = z['source_force'][index]
        self.saved_excitation[:] = z['excitation'][index]
        self.reset()

    def reset(self):
        m, d, b = self.model, self.data, self.backend
        m.actuator_biasprm[:] = self.saved_bias
        mujoco.mj_setState(m, d, self.saved_state, self.mask)
        b.bank.state[:] = self.saved_fiber
        b.bank.force[:] = self.saved_force
        b.bank.bounds[:] = 0
        b.knee_max[:] = d.qpos[b.knees]
        self.excitation[:] = self.saved_excitation
        mujoco.mj_forward(m, d)

    def commit_snapshot(self):
        """Carry an accepted physical boundary, including every latent state."""
        mujoco.mj_getState(self.model,self.data,self.saved_state,self.mask)
        self.saved_bias[:] = self.model.actuator_biasprm
        self.saved_fiber[:] = self.backend.bank.state
        self.saved_force[:] = self.backend.bank.force
        self.saved_excitation[:] = self.data.ctrl

    def rollout(self, command, steps, recorded=None, observer=None):
        if steps > len(self.activation_distance):
            raise ValueError('Diagnostic window exceeds preallocated branch buffers')
        self.reset()
        for k in range(steps):
            if recorded is None:
                self.excitation += self.beta*(command-self.excitation)
                self.data.ctrl[:] = self.excitation
            else:
                self.data.ctrl[:] = recorded[k]
            self.activation_distance[k] = self.data.ctrl-self.data.act
            self.backend.step()
            self.fiber_velocity[k] = (self.backend.bank.state-self.backend.old_state)/self.model.opt.timestep
            self.physics_steps += 1
            if self.source_rom is not None:
                address, low, high = self.source_rom
                if np.any(self.data.qpos[address]<low-1e-10) or np.any(self.data.qpos[address]>high+1e-10):
                    raise RuntimeError('Frozen source ROM in augmented rollout')
            if observer is not None:
                observer()
        d = self.data
        return np.r_[d.qpos, d.qvel, d.act, self.backend.bank.state, d.ctrl]


class NativeContactKinematics:
    """Relative material-point slip and the actual native contact force cone."""
    def __init__(self,model,data):
        self.model,self.data = model,data
        self.ownership = {model.body(n+'_'+s).id:i for i,s in enumerate(('r','l')) for n in ('calcn','toes')}
        self.force = np.zeros(6); self.jp = np.zeros((3,model.nv)); self.jr = np.zeros_like(self.jp)
        self.velocity = np.zeros((2,3)); self.angular = np.zeros_like(self.velocity)
    def measure(self,branches=False):
        m,d = self.model,self.data; normals = np.zeros(2); slip = np.zeros(2)
        body_normals=np.zeros(4);body_counts=np.zeros(4,dtype=int);signature=[]
        minimum_normal = np.inf; maximum_violation = 0.; maximum_utilization = 0.; contact_count = 0
        for k in range(d.ncon):
            c = d.contact[k]; bodies = (int(m.geom_bodyid[c.geom1]),int(m.geom_bodyid[c.geom2]))
            sides = {self.ownership[b] for b in bodies if b in self.ownership}
            if not sides: continue
            if len(sides)!=1: raise RuntimeError('cross_foot_contact')
            if c.efc_address<0: continue
            side = sides.pop(); mujoco.mj_contactForce(m,d,k,self.force)
            fn = float(self.force[0]); minimum_normal = min(minimum_normal,fn); normals[side] += fn; contact_count += 1
            if branches:
                body=next(body for body in bodies if body in self.ownership);part=2*side+int(m.body(body).name.startswith('toes'))
                body_normals[part]+=fn;body_counts[part]+=1;signature.append((part,int(c.geom1),int(c.geom2),int(c.dim)))
            for j,body in enumerate(bodies):
                mujoco.mj_jac(m,d,self.jp,self.jr,c.pos,body)
                self.velocity[j] = self.jp@d.qvel; self.angular[j] = self.jr@d.qvel
            relative = c.frame.reshape(3,3)@(self.velocity[1]-self.velocity[0])
            if fn>.001*float(m.body_mass.sum()*9.81):
                slip[side] = max(slip[side],float(np.linalg.norm(relative[1:])))
            dimension = int(c.dim)
            scaled = self.force[1:dimension]/c.friction[:dimension-1]
            tangent_equivalent = float(np.sum(abs(scaled)) if m.opt.cone==mujoco.mjtCone.mjCONE_PYRAMIDAL else np.linalg.norm(scaled))
            maximum_violation = max(maximum_violation,tangent_equivalent-fn)
            if fn>.001*float(m.body_mass.sum()*9.81):
                maximum_utilization = max(maximum_utilization,tangent_equivalent/fn)
        result=dict(normal_N=normals,slip_mps=slip,minimum_normal_N=float(minimum_normal) if contact_count else 0.,
                    maximum_friction_cone_violation_N=maximum_violation,maximum_loaded_cone_utilization=maximum_utilization,
                    contacts=contact_count)
        if branches:
            import zlib
            result.update(normal_by_body_N=body_normals,contact_count_by_body=body_counts,signature_u32=zlib.crc32(repr(sorted(signature)).encode()))
        return result


def validate_contact_kinematics(model_folder,trace_folders,out):
    """Independent contact material-point FK and contact virtual-work checks."""
    m = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb')); d = mujoco.MjData(m)
    plus,minus = mujoco.MjData(m),mujoco.MjData(m)
    jacp = np.zeros((3,m.nv)); jacr = np.zeros_like(jacp); f = np.zeros(6); tau = np.zeros(m.nv)
    m.actuator_gaintype[:] = mujoco.mjtGain.mjGAIN_FIXED; m.actuator_gainprm[:] = 0.
    m.actuator_biastype[:] = mujoco.mjtBias.mjBIAS_AFFINE; m.actuator_biasprm[:] = 0.
    maximum_velocity_error = maximum_work_error = maximum_cone = maximum_projection_error = 0.; minimum_normal = np.inf; checks = 0
    if m.nq!=m.nv: raise ValueError('Scalar captured model required for independent FK differences')
    fixture_foot_angles = []
    for folder in trace_folders:
        z = np.load(folder/'final_continuous_trace.npz')
        for index in np.unique(np.linspace(0,len(z['t'])-1,5).astype(int)):
            d.qpos[:] = z['qpos'][index]; d.qvel[:] = z['qvel'][index]; d.act[:] = z['activation'][index]
            d.xfrc_applied[m.body('pelvis').id,2] = .6*float(m.body_mass.sum()*9.81)
            m.actuator_biasprm[:,0] = -z['source_force_N'][index]
            mujoco.mj_forward(m,d)
            probe = NativeContactKinematics(m,d); values = probe.measure()
            minimum_normal = min(minimum_normal,values['minimum_normal_N']); maximum_cone = max(maximum_cone,values['maximum_friction_cone_violation_N'])
            h = 1e-4
            plus.qpos[:] = d.qpos+h*d.qvel; minus.qpos[:] = d.qpos-h*d.qvel
            mujoco.mj_fwdPosition(m,plus); mujoco.mj_fwdPosition(m,minus)
            for k in range(d.ncon):
                c = d.contact[k]; bodies = (int(m.geom_bodyid[c.geom1]),int(m.geom_bodyid[c.geom2]))
                if not any(b in probe.ownership for b in bodies) or c.efc_address<0: continue
                relative = np.zeros(3); angular = np.zeros(3)
                for sign,body in zip((-1.,1.),bodies):
                    local = d.xmat[body].reshape(3,3).T@(c.pos-d.xpos[body])
                    measured = ((plus.xpos[body]+plus.xmat[body].reshape(3,3)@local)-(minus.xpos[body]+minus.xmat[body].reshape(3,3)@local))/(2*h)
                    mujoco.mj_jac(m,d,jacp,jacr,c.pos,body)
                    maximum_velocity_error = max(maximum_velocity_error,float(np.max(abs(measured-jacp@d.qvel))))
                    relative += sign*(jacp@d.qvel); angular += sign*(jacr@d.qvel)
                mujoco.mj_contactForce(m,d,k,f); rotation = c.frame.reshape(3,3)
                world_force = rotation.T@f[:3]; world_torque = rotation.T@f[3:]
                tau[:] = 0.
                for sign,body in zip((-1.,1.),bodies):
                    mujoco.mj_applyFT(m,d,sign*world_force,sign*world_torque,c.pos,body,tau)
                contact_force = np.zeros(d.nefc); native_projection = np.zeros(m.nv)
                count = 1 if c.dim==1 else 2*(int(c.dim)-1) if m.opt.cone==mujoco.mjtCone.mjCONE_PYRAMIDAL else int(c.dim)
                contact_force[c.efc_address:c.efc_address+count] = d.efc_force[c.efc_address:c.efc_address+count]
                mujoco.mj_mulJacTVec(m,d,native_projection,contact_force)
                maximum_projection_error = max(maximum_projection_error,float(np.max(abs(tau-native_projection))))
                power = float(world_force@relative+world_torque@angular)
                maximum_work_error = max(maximum_work_error,abs(power-float(tau@d.qvel)))
                checks += 1
            fixture_foot_angles.append([float(d.qpos[m.joint(n+'_'+s).qposadr[0]]) for s in ('r','l') for n in ('subtalar_angle','mtp_angle')])
    result = dict(material_point_FK_velocity_error_mps=maximum_velocity_error,contact_virtual_work_error_W=maximum_work_error,
                  native_contact_generalized_force_projection_error=maximum_projection_error,
                  minimum_native_normal_N=minimum_normal,maximum_native_friction_cone_violation_N=maximum_cone,contact_checks=checks,
                  foot_teacher_posture_degrees=np.degrees(fixture_foot_angles).tolist(),
                  reference='https://mujoco.readthedocs.io/en/2.3.7/computation.html#contact',
                  scope='Read-only independent kinematics/force-coordinate validation on saved loaded states; not a new physical trajectory.')
    result['passed'] = bool(checks>0 and maximum_velocity_error<1e-7 and maximum_work_error<1e-8 and maximum_projection_error<1e-7 and minimum_normal>=-1e-8 and maximum_cone<1e-6)
    out.mkdir(parents=True,exist_ok=True); target = out/'CONTACT_KINEMATICS_VALIDATION.json'
    if target.exists(): raise ValueError('Preserve contact validation')
    target.write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))
    if not result['passed']: raise SystemExit(1)


def stationary_group_inverse(model_folder, start_folder, end_folder, spec_path, out, verify_folder=None, physics_dt=None, continuation_seed=None, verify_continuation=None,rolling_stance=False,extend_transfer_from=None):
    """Bounded nonlinear load-transfer test from an accepted native equilibrium.

    Group amplitudes multiply saved endpoint-state differences, not independently
    learned leg gains. Every residual evaluation advances full native contact,
    activation, filtering, and source fibers. No force Jacobian is substituted.
    """
    from time import perf_counter
    from scipy.optimize import minimize, brentq
    from scipy.linalg import expm
    from thelen_dynamic import ConstitutiveBackend, balance_and_slope, force_velocity
    from source_force_comparison import load_recruitment
    if (out/'STATIONARY_GROUP_INVERSE.json').exists():
        raise ValueError('Preserve existing nonlinear inverse result')
    start, end = (load_recruitment(f) for f in (start_folder, end_folder))
    reports = [json.loads((f/'SOURCE_NATIVE_STANCE.json').read_text()) for f in (start_folder, end_folder)]
    if not all(r.get('accepted_stationary_endpoint', False) for r in reports):
        raise ValueError('Both source endpoints must pass existing stationarity gates')
    endpoint_normal = np.asarray(reports[1]['metrics']['last_half_normal_min_BW'])
    if endpoint_normal.shape != (2,) or np.min(endpoint_normal)>.025 or np.max(endpoint_normal)<.05:
        raise ValueError('An accepted single-support endpoint is required')
    if continuation_seed is None and (endpoint_normal[0]>.025 or endpoint_normal[1]<.05):
        raise ValueError('This declared right-to-left task requires an accepted right-unloaded endpoint')
    config = json.loads((start_folder/'MANIFEST.json').read_text())
    support = float(config['support_fraction'])
    for key in ('parameters', 'muscles', 'force_channel_enabled'):
        if not np.array_equal(start[key], end[key]):
            raise ValueError('Endpoint source parameters or actuation policies differ')
    if not np.all(start['force_channel_enabled']):
        raise ValueError('This diagnostic explicitly uses the qualified source92 endpoints')
    m = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
    if physics_dt is not None:
        m.opt.timestep = physics_dt
    d = mujoco.MjData(m); p = start['parameters']; command = start['activation'].copy()
    names = start['muscles'].tolist(); dt = float(m.opt.timestep); bw = float(m.body_mass.sum()*9.81)
    d.qpos[:] = start['qpos']; d.act[:] = command; d.ctrl[:] = command
    d.xfrc_applied[m.body('pelvis').id, 2] = support*bw
    mujoco.mj_forward(m, d)
    source = ConstitutiveBackend(m, d, p)
    m.actuator_biasprm[:, 0] = -source.bank.force
    mujoco.mj_forward(m, d)
    initial_force = source.bank.force.copy()
    initial_acceleration = d.qacc.copy()
    source_force_error = float(np.max(abs(initial_force-start['source_force'])))
    if source_force_error > 1e-8:
        raise ValueError('Qualified double-support source force was not reproduced')
    families = dict(hip_flexor_mono=('psoas','iliacus'), hip_extensor_mono=('glut_max1','glut_max2','glut_max3'),
                    hip_abductor_mono=tuple(x+str(k) for x in ('glut_med','glut_min') for k in (1,2,3)),
                    hip_adductor_mono=('add_long','add_brev','add_mag1','add_mag2','add_mag3','pect'),
                    knee_extensor_mono=('vas_lat','vas_med','vas_int'), knee_flexor_mono=('bifemsh',),
                    biarticular_HAM=('semimem','semiten','bifemlh'), biarticular_RF=('rect_fem',),
                    biarticular_GAS=('med_gas','lat_gas'), ankle_plantarflexor_mono=('soleus',),
                    ankle_dorsiflexor_mono=('tib_ant',))
    delta = end['activation']-command; columns = []; group_names = []; limits = []; global_limits = []
    for label, members in families.items():
        column = np.array([delta[j] if n[:-2] in members else 0. for j,n in enumerate(names)])
        maximum = float(np.max(abs(column)))
        if maximum > 1e-12:
            columns.append(column/maximum); group_names.append(label); limits.append(min(.02, maximum))
            global_limits.append(maximum)
    B = np.array(columns).T; limits = np.array(limits)
    # Derive a fixed response time from the native activation law and the source
    # constitutive velocity tangent, retaining the engineering 50 ms filter.
    pole_rows = []
    selected = np.flatnonzero(np.any(B != 0., axis=1))
    for j in selected:
        a = command[j]; scale = .5+1.5*a
        tau_a = float(m.actuator_dynprm[j, 0]*scale if delta[j]>0. else m.actuator_dynprm[j, 1]/scale)
        s = source.bank.state[j]; fiber = np.hypot(s, p[j,3]); co = s/fiber; nl = fiber/p[j,1]
        _, stiffness, _, _, _ = balance_and_slope(s,s,d.actuator_length[j],a,p[j],dt,False)
        fl = np.exp(-(nl-1.)**2/.5)
        dfv = max(force_velocity(-1e-14,a,p[j,6],p[j,7],p[j,8])[1], force_velocity(1e-14,a,p[j,6],p[j,7],p[j,8])[1])
        tau_f = max(1e-9, (max(.01,a)*fl*dfv+p[j,5])*co*co/(p[j,1]*stiffness))
        poles = np.array([.05, tau_a, tau_f])
        # Matrix exponential handles equal poles without unstable subtraction.
        A = np.diag(-1./poles); A[1,0] = 1./poles[1]; A[2,1] = 1./poles[2]
        def response(t):
            return float(1.-(expm(A*t)@np.ones(3))[2])
        t90 = brentq(lambda t: response(t)-.9, 0., 20.*poles.sum())
        pole_rows.append(dict(muscle=names[j], activation_tau_s=tau_a, fiber_tau_s=tau_f,
                              filter_tau_s=.05, force_response_fraction_10ms=response(.01), t90_s=t90))
    derived_horizon = np.ceil(max(v['t90_s'] for v in pole_rows)/dt)*dt
    if derived_horizon > .5:
        raise ValueError('Derived horizon exceeds the declared bounded half-second diagnostic')
    allocated_horizon = max(.302,derived_horizon) if continuation_seed is not None else derived_horizon
    runner = AugmentedTransition(m, p, .05, max_steps=max(16, int(round(allocated_horizon/dt))))
    # Construction changes force-law arrays on the cloned model; restore the
    # exact already validated source-force state before snapshotting it.
    rd, rm = runner.data, runner.model
    rd.qpos[:] = d.qpos; rd.qvel[:] = d.qvel; rd.act[:] = d.act; rd.ctrl[:] = d.ctrl
    rd.xfrc_applied[:] = d.xfrc_applied; rm.actuator_biasprm[:] = m.actuator_biasprm
    mujoco.mj_forward(rm, rd)
    mujoco.mj_getState(rm, rd, runner.saved_state, runner.mask)
    runner.saved_bias[:] = rm.actuator_biasprm
    runner.saved_fiber[:] = source.bank.state; runner.saved_force[:] = initial_force
    runner.saved_excitation[:] = command
    spec = json.loads(spec_path.read_text())
    address = np.array([v['qpos_address'] for v in spec['joints'].values()])
    low = np.array([v['native_bounds_rad'][0] for v in spec['joints'].values()])
    high = np.array([v['native_bounds_rad'][1] for v in spec['joints'].values()])
    runner.source_rom = (address, low, high)
    feet = [rm.body('calcn_'+s).id for s in ('r','l')]
    ownership = {rm.body(n+'_'+s).id:i for i,s in enumerate(('r','l')) for n in ('calcn','toes')}
    force6 = np.zeros(6); jp = np.zeros((3,rm.nv)); jr = np.zeros_like(jp)
    def measurements():
        normals = np.zeros(2); slip = np.zeros(2)
        for k in range(rd.ncon):
            contact = rd.contact[k]
            sides = {ownership[b] for b in (int(rm.geom_bodyid[contact.geom1]),int(rm.geom_bodyid[contact.geom2])) if b in ownership}
            if not sides:
                continue
            if len(sides) != 1:
                raise ValueError('Unexpected cross-foot contact in supported-stance diagnostic')
            side = sides.pop(); mujoco.mj_contactForce(rm,rd,k,force6)
            normal = max(0.,float(force6[0])); normals[side] += normal
            if normal > .001*bw:
                body = next(b for b in (int(rm.geom_bodyid[contact.geom1]),int(rm.geom_bodyid[contact.geom2])) if b in ownership)
                mujoco.mj_jac(rm,rd,jp,jr,contact.pos,body)
                velocity = jp@rd.qvel
                slip[side] = max(slip[side],float(np.linalg.norm(velocity[:2])))
        return normals, slip
    runner.reset(); initial_normal, _ = measurements()
    # One percent BW transfer is an explicit diagnostic task, not a replacement
    # for the existing complete unloading, clearance, or touchdown gates.
    target_normal = initial_normal+np.array([-.01,.01])*bw
    if np.min(target_normal)<0.:
        raise ValueError('Diagnostic expects loaded double support')
    q0 = rd.qpos.copy()
    if continuation_seed is not None:
        continue_group_transfer(runner,command,B,group_names,np.asarray(global_limits),q0,d.xpos[feet].copy(),
                                feet,ownership,measurements,initial_normal,endpoint_normal,spec,address,
                                continuation_seed,out,support,start_folder,end_folder,verify_continuation,rolling_stance,extend_transfer_from)
        return
    if verify_folder is not None:
        accepted = json.loads((verify_folder/'STATIONARY_GROUP_INVERSE.json').read_text())
        prior = np.load(verify_folder/'bounded_states.npz')
        if accepted['group_names'] != group_names or not np.array_equal(prior['group_basis'], B):
            raise ValueError('Saved functional basis changed')
        ref = np.load(spec['source_target'])
        speed_limit = np.maximum(np.max(abs(ref['qvel'][:,address]),axis=0),1e-8)
        acceleration_limit = np.maximum(np.max(abs(ref['qacc'][:,address]),axis=0),1e-8)
        geometry = []
        for g in range(rm.ngeom):
            if int(rm.geom_bodyid[g]) in ownership and rm.geom_type[g] == mujoco.mjtGeom.mjGEOM_MESH:
                mesh = int(rm.geom_dataid[g]); first = rm.mesh_vertadr[mesh]; count = rm.mesh_vertnum[mesh]
                geometry.append((ownership[int(rm.geom_bodyid[g])],g,rm.mesh_vert[first:first+count].copy()))
        pg = rm.geom('ground-plane').id
        results = []
        for index,row in enumerate(accepted['rows']):
            horizon = float(row['horizon_s']); steps = round(horizon/dt)
            if abs(steps*dt-horizon)>1e-12:
                raise ValueError('Saved horizon is not aligned to the declared timestep')
            maximum_slip = np.zeros(2); minimum_gap = np.full(2,np.inf)
            maximum_root = np.zeros(2); maximum_rates = np.zeros(2)
            def observe():
                _, slip = measurements()
                maximum_slip[:] = np.maximum(maximum_slip,slip)
                origin = rd.geom_xpos[pg]; normal = rd.geom_xmat[pg].reshape(3,3)[:,2]
                for side,g,vertices in geometry:
                    gap = (vertices@rd.geom_xmat[g].reshape(3,3).T+rd.geom_xpos[g]-origin)@normal
                    minimum_gap[side] = min(minimum_gap[side],float(gap.min()))
                maximum_root[:] = np.maximum(maximum_root,[np.max(abs(rd.qpos[:3]-q0[:3])),np.rad2deg(np.max(abs(rd.qpos[3:6]-q0[3:6])))])
                maximum_rates[:] = np.maximum(maximum_rates,[np.max(abs(rd.qvel[address])/speed_limit),np.max(abs(rd.qacc[address])/acceleration_limit)])
            u = command+B@np.array(row['best_group_correction'])
            runner.rollout(u,steps,observer=observe)
            normal, _ = measurements()
            error = {key:float(np.max(abs(value-prior[str(index)+'_'+key]))) for key,value in dict(qpos=rd.qpos,qvel=rd.qvel,activation=rd.act,fiber_state=runner.backend.bank.state).items()}
            result = dict(horizon_s=horizon,timestep_s=dt,final_state_difference_from_original=error,
                          normal_load_BW=(normal/bw).tolist(), normal_error_BW=((normal-target_normal)/bw).tolist(),
                          maximum_loaded_slip_mps=maximum_slip.tolist(),minimum_visual_mesh_gap_m=minimum_gap.tolist(),
                          maximum_root_translation_m=float(maximum_root[0]),maximum_root_rotation_deg=float(maximum_root[1]),
                          maximum_source_speed_ratio=float(maximum_rates[0]),maximum_source_acceleration_ratio=float(maximum_rates[1]),
                          maximum_knee_extension_deg=np.rad2deg(runner.backend.knee_max).tolist(),
                          backend=runner.backend.diagnostics())
            result['window_physics_gates_pass'] = bool(max(maximum_slip)<=.005 and min(minimum_gap)>=-.002 and maximum_root[0]<=.005 and maximum_root[1]<=2. and max(maximum_rates)<=1. and max(runner.backend.knee_max)<=0.)
            result['load_transfer_task_pass'] = bool(np.max(abs(normal-target_normal))/bw<=.002)
            results.append(result)
        out.mkdir(parents=True,exist_ok=True)
        target = out/'STATIONARY_GROUP_REPLAY.json'
        if target.exists():
            raise ValueError('Preserve existing stationary-group replay')
        report = dict(source=str(verify_folder),timestep_s=dt,rows=results,
                      scope='Fresh saved-command replay,all-window physical checks; no fitting or physical changes. This1%BW load-transfer diagnostic is not a complete unloading/clearance teacher.')
        target.write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2)); return
    rows = []; saved = []; started = perf_counter()
    for horizon in (.01, float(derived_horizon)):
        if rows and rows[-1]['task_feasible']:
            break
        steps = int(round(horizon/dt)); history = []; best = [np.inf, None, None]
        def evaluate(x, store=True):
            if np.any(x < -1e-12) or np.any(x > limits+1e-12):
                return 1e12
            u = command+B@x
            if np.any(u<.01-1e-12) or np.any(u>.98+1e-12):
                return 1e12
            try:
                runner.rollout(u,steps)
                normal, slip = measurements()
                foot_error = rd.xpos[feet]-d.xpos[feet]
                values = dict(normal_load_BW=(normal/bw).tolist(), target_normal_load_BW=(target_normal/bw).tolist(),
                              normal_error_BW=((normal-target_normal)/bw).tolist(), loaded_slip_mps=slip.tolist(),
                              foot_position_error_m=float(np.max(abs(foot_error))), root_translation_error_m=float(np.max(abs(rd.qpos[:3]-q0[:3]))),
                              root_rotation_error_deg=float(np.rad2deg(np.max(abs(rd.qpos[3:6]-q0[3:6])))),
                              knee_max_deg=float(np.rad2deg(np.max(runner.backend.knee_max))))
                components = dict(normal_load=float(np.mean(((normal-target_normal)/(.01*bw))**2)),
                                  foot_position=float(np.mean((foot_error/.0002)**2)),
                                  loaded_slip=float(np.mean((slip/.005)**2)),
                                  root_translation=float(np.mean(((rd.qpos[:3]-q0[:3])/.005)**2)),
                                  root_rotation=float(np.mean(((rd.qpos[3:6]-q0[3:6])/np.deg2rad(2.))**2)))
                cost = float(sum(components.values()))
                values['normalized_component_mean_squares'] = components
                values['guard'] = ''
            except RuntimeError as error:
                cost = 1e10; values = dict(guard=str(error))
            if store:
                history.append(dict(x=x.tolist(),cost=cost,**values))
                if cost<best[0]:
                    best[:] = [cost,x.copy(),values]
            return cost
        zero = np.zeros(len(limits)); evaluate(zero); baseline = history[-1].copy()
        for j in range(len(limits)):
            probe = zero.copy(); probe[j] = limits[j]; evaluate(probe)
        # Native branch changes are evaluated directly. No smooth finite-
        # difference Jacobian or independently fitted leg gains are assumed.
        seed = best[1].copy()
        constraints = [dict(type='ineq',fun=lambda x:x),dict(type='ineq',fun=lambda x:limits-x)]
        fit = minimize(evaluate,seed,method='COBYLA',constraints=constraints,
                       options=dict(rhobeg=.005,tol=1e-5,maxiter=100,catol=1e-10))
        x = best[1]; evaluate(x,store=False); final = best[2]
        feasible = not final['guard'] and max(abs(np.array(final['normal_error_BW'])))<=.002 and max(final['loaded_slip_mps'])<=.005 and final['foot_position_error_m']<=.0002 and final['root_translation_error_m']<=.005 and final['root_rotation_error_deg']<=2. and final['knee_max_deg']<=0.
        rows.append(dict(horizon_s=horizon,baseline=baseline,best=final,best_group_correction=x.tolist(),
                         cost_reduction_fraction=float(1.-best[0]/baseline['cost']), evaluations=len(history),task_feasible=bool(feasible),search_history=history))
        saved.append(dict(horizon_s=horizon,command=command+B@x,qpos=rd.qpos.copy(),qvel=rd.qvel.copy(),activation=rd.act.copy(),fiber_state=runner.backend.bank.state.copy()))
    result = dict(start=str(start_folder),target=str(end_folder),source92=True,support_fraction=support,
                  source_command_sha256=hashlib.sha256(command.tobytes()).hexdigest(),initial_native_acceleration_max=float(np.max(abs(initial_acceleration))),
                  initial_source_force_error_N=source_force_error,
                  group_names=group_names,group_upper_bounds=limits.tolist(),group_rule='One common amplitude per functional family multiplies its saved bilateral endpoint-state difference. Other muscles remain at qualified double-support commands; RF,HAM,GAS are separate; SAR/GRAC unchanged.',
                  task='Transfer1%BW from right to left while retaining loaded stationary foot geometry. Not a swing/teacher acceptance.',
                  component_scales=dict(normal_BW=.01,foot_position_m=.0002,loaded_contact_slip_mps=.005,root_translation_m=.005,root_rotation_deg=2.),
                  analytic_response=dict(method='Three-pole filter/native-activation/source-velocity-tangent small-signal fixed-MTU response; larger one-sided FV slope. Excludes changing mechanical geometry/contact,not a whole-model settling-time prediction.',derived_next_horizon_s=float(derived_horizon),muscles=pole_rows),
                  rows=rows,native_physics_steps=runner.physics_steps,wall_s=perf_counter()-started,
                  scope='Offline nonlinear augmented command correction only. No source physics,endpoint,gain,contact,force-floor change; no live controller or Hunt C fit.')
    out.mkdir(parents=True,exist_ok=True)
    (out/'STATIONARY_GROUP_INVERSE.json').write_text(json.dumps(result,indent=2))
    np.savez_compressed(out/'bounded_states.npz',group_basis=B,muscles=names,
                        **{str(k)+'_'+key:value for k,row in enumerate(saved) for key,value in row.items()})
    print(json.dumps({k:v for k,v in result.items() if k not in ('rows','analytic_response')},indent=2))
    print(json.dumps([{k:v for k,v in row.items() if k!='search_history'} for row in rows],indent=2))


def continue_group_transfer(runner,origin,B,groups,global_limits,q0,foot0,feet,ownership,measurements,
                            initial_normal,endpoint_normal,spec,address,seed_folder,out,support,start_folder,end_folder,replay_folder=None,rolling_stance=False,extend_transfer_from=None):
    """Continuous1/3/5%BW sequence; each accepted boundary carries all states."""
    from time import perf_counter
    from scipy.optimize import minimize
    output_name = 'CONTINUOUS_GROUP_REPLAY.json' if replay_folder is not None else 'CONTINUOUS_GROUP_TRANSFER.json'
    if (out/output_name).exists():
        raise ValueError('Preserve existing continuous transfer')
    seed = json.loads((seed_folder/'STATIONARY_GROUP_INVERSE.json').read_text())
    if seed['group_names'] != groups:
        raise ValueError('Functional families differ from preserved1% reference')
    accepted_seed = next(row for row in seed['rows'] if row['task_feasible'])
    horizon = float(accepted_seed['horizon_s'])
    if abs(horizon-.302)>1e-12:
        raise ValueError('Continuation is restricted to the verified302ms horizon')
    m,d = runner.model,runner.data; dt = float(m.opt.timestep); steps = round(horizon/dt)
    if abs(steps*dt-horizon)>1e-12:
        raise ValueError('Fixed continuation horizon is not timestep aligned')
    bw = float(m.body_mass.sum()*9.81); pg = m.geom('ground-plane').id
    unloading = int(np.argmin(endpoint_normal)); direction = np.ones(2); direction[unloading] = -1.
    ref = np.load(spec['source_target'])
    speed_limit = np.maximum(np.max(abs(ref['qvel'][:,address]),axis=0),1e-8)
    acceleration_limit = np.maximum(np.max(abs(ref['qacc'][:,address]),axis=0),1e-8)
    geometry = []
    for g in range(m.ngeom):
        if int(m.geom_bodyid[g]) in ownership and m.geom_type[g]==mujoco.mjtGeom.mjGEOM_MESH:
            mesh = int(m.geom_dataid[g]); first = m.mesh_vertadr[mesh]; count = m.mesh_vertnum[mesh]
            geometry.append((ownership[int(m.geom_bodyid[g])],g,m.mesh_vert[first:first+count].copy()))
    contact_probe = NativeContactKinematics(m,d) if rolling_stance else None
    foot_posture_names = [n+'_'+s for s in ('r','l') for n in ('subtalar_angle','mtp_angle')]
    foot_posture_address = np.array([m.joint(n).qposadr[0] for n in foot_posture_names])
    # This is the pre-existing full-native static teacher envelope, not the
    # original OpenSim mathematical +/-90 degree unclamped coordinate range.
    foot_posture_limit = np.deg2rad(10.)
    trace = []; stages = []; previous = np.zeros(len(groups)); previous_fraction = 0.
    replay = None
    extension = None; extension_trace = None
    fractions = (.01,.03,.05)
    if extend_transfer_from is not None:
        if not rolling_stance:
            raise ValueError('Extended transfer is a distinct rolling-stance task only')
        extension = json.loads((extend_transfer_from/'CONTINUOUS_GROUP_TRANSFER.json').read_text())
        extension_trace = np.load(extend_transfer_from/'final_continuous_trace.npz')
        if not extension.get('rolling_stance',False) or extension['groups']!=groups or extension['unloading_side']!=('r','l')[unloading]:
            raise ValueError('Extension changes the accepted phase or functional policy')
        if not np.array_equal(extension_trace['group_basis'],B) or not np.array_equal(extension_trace['parameters'],runner.backend.parameters):
            raise ValueError('Extension changes the accepted basis or muscle parameters')
        old_fractions = [row['fraction_BW'] for row in extension['stages'] if row['passed']]
        if old_fractions==[.01,.03,.05]:
            fractions = tuple(old_fractions)+(.075,.10)
        elif old_fractions==[.01,.03,.05,.075,.10]:
            fractions = tuple(old_fractions)+(.125,.15)
        else:
            raise ValueError('Extension requires an explicitly supported accepted5% or10% prefix')
    if replay_folder is not None:
        replay = json.loads((replay_folder/'CONTINUOUS_GROUP_TRANSFER.json').read_text())
        if replay['groups']!=groups or replay['unloading_side']!=('r','l')[unloading]:
            raise ValueError('Saved continuous-transfer phase or groups changed')
        if bool(replay.get('rolling_stance',False))!=rolling_stance:
            raise ValueError('Saved stance task variant differs')
        fractions = tuple(row['fraction_BW'] for row in replay['stages'] if row['passed'])
        if not fractions:
            raise ValueError('No accepted continuous prefix to replay')
    started = perf_counter()
    for fraction in fractions:
        low = np.maximum(0.,previous-.02); high = np.minimum(global_limits,previous+.02)
        target = initial_normal+direction*fraction*bw
        target_start = initial_normal+direction*previous_fraction*bw
        if min(target)<0.:
            raise ValueError('Declared transfer exceeds initial loaded-foot force')
        history = []; best = [np.inf,None,None]; current_trace = []
        def evaluate(x,store=True,record=False):
            if np.any(x<low-1e-12) or np.any(x>high+1e-12):
                return 1e12
            u = origin+B@x
            if min(u)<.01-1e-12 or max(u)>.98+1e-12:
                return 1e12
            maximum_slip = np.zeros(2); minimum_gap = np.full(2,np.inf)
            maxima = dict(foot_position_m=0.,root_translation_m=0.,root_rotation_deg=0.,source_speed_ratio=0.,
                          source_acceleration_ratio=0.,fiber_state_speed_per_l0_s=0.,fiber_length_speed_per_l0_s=0.)
            contact_maxima = dict(minimum_normal_N=np.inf,friction_cone_violation_N=0.,loaded_cone_utilization=0.,absolute_foot_posture_deg=0.)
            path_cost = []; sample = []; start_time = float(d.time)
            # The solver may have left a trial at its end: use the immutable
            # committed native boundary time,not that disposable current time.
            runner.reset(); start_time = float(d.time)
            def observe():
                contact = contact_probe.measure() if rolling_stance else None
                normal,slip = (contact['normal_N'],contact['slip_mps']) if rolling_stance else measurements()
                maximum_slip[:] = np.maximum(maximum_slip,slip)
                if rolling_stance:
                    contact_maxima['minimum_normal_N'] = min(contact_maxima['minimum_normal_N'],contact['minimum_normal_N'])
                    contact_maxima['friction_cone_violation_N'] = max(contact_maxima['friction_cone_violation_N'],contact['maximum_friction_cone_violation_N'])
                    contact_maxima['loaded_cone_utilization'] = max(contact_maxima['loaded_cone_utilization'],contact['maximum_loaded_cone_utilization'])
                    contact_maxima['absolute_foot_posture_deg'] = max(contact_maxima['absolute_foot_posture_deg'],float(np.rad2deg(np.max(abs(d.qpos[foot_posture_address])))))
                ground_normal = d.geom_xmat[pg].reshape(3,3)[:,2]; ground_origin = d.geom_xpos[pg]
                instantaneous_gap = np.full(2,np.inf)
                for side,g,vertices in geometry:
                    local_normal = d.geom_xmat[g].reshape(3,3).T@ground_normal
                    gap = np.min(vertices@local_normal)+(d.geom_xpos[g]-ground_origin)@ground_normal
                    minimum_gap[side] = min(minimum_gap[side],float(gap))
                    instantaneous_gap[side] = min(instantaneous_gap[side],float(gap))
                maxima['foot_position_m'] = max(maxima['foot_position_m'],float(np.max(abs(d.xpos[feet]-foot0))))
                maxima['root_translation_m'] = max(maxima['root_translation_m'],float(np.max(abs(d.qpos[:3]-q0[:3]))))
                maxima['root_rotation_deg'] = max(maxima['root_rotation_deg'],float(np.rad2deg(np.max(abs(d.qpos[3:6]-q0[3:6])))))
                maxima['source_speed_ratio'] = max(maxima['source_speed_ratio'],float(np.max(abs(d.qvel[address])/speed_limit)))
                maxima['source_acceleration_ratio'] = max(maxima['source_acceleration_ratio'],float(np.max(abs(d.qacc[address])/acceleration_limit)))
                fiber_state_speed = (runner.backend.bank.state-runner.backend.old_state)/(dt*runner.backend.parameters[:,1])
                pennation_cosine = runner.backend.bank.state/np.hypot(runner.backend.bank.state,runner.backend.parameters[:,3])
                maxima['fiber_state_speed_per_l0_s'] = max(maxima['fiber_state_speed_per_l0_s'],float(np.max(abs(fiber_state_speed))))
                maxima['fiber_length_speed_per_l0_s'] = max(maxima['fiber_length_speed_per_l0_s'],float(np.max(abs(fiber_state_speed*pennation_cosine))))
                phase = np.clip((float(d.time)-start_time)/horizon,0.,1.)
                blend = phase**3*(10.+phase*(-15.+6.*phase))
                desired = target_start+blend*(target-target_start)
                path_cost.append(float(np.mean(((normal-desired)/(.01*bw))**2)))
                if record:
                    sample.append(dict(t=float(d.time),qpos=d.qpos.copy(),qvel=d.qvel.copy(),qacc=d.qacc.copy(),
                                       activation=d.act.copy(),fiber_state=runner.backend.bank.state.copy(),
                                       excitation=d.ctrl.copy(),command=u.copy(),source_force_N=runner.backend.bank.force.copy(),normal_load_BW=normal/bw,
                                       mesh_gap_m=instantaneous_gap.copy()))
                reasons = []
                if max(maximum_slip)>.005: reasons.append('loaded_contact_slip')
                if min(minimum_gap)<-.002: reasons.append('visual_mesh_penetration')
                if not rolling_stance and maxima['foot_position_m']>.0002: reasons.append('absolute_loaded_foot_position')
                if rolling_stance:
                    if contact_maxima['minimum_normal_N'] < -1e-8: reasons.append('unilateral_normal_force')
                    if contact_maxima['friction_cone_violation_N'] > 1e-6: reasons.append('native_friction_cone')
                    if np.max(abs(d.qpos[foot_posture_address]))>foot_posture_limit: reasons.append('source_teacher_foot_posture')
                if maxima['root_translation_m']>.005 or maxima['root_rotation_deg']>2.: reasons.append('absolute_root_pose')
                if maxima['source_speed_ratio']>1. or maxima['source_acceleration_ratio']>1.: reasons.append('source_joint_rate')
                if min(d.act)<.01-1e-12 or max(d.act)>.98+1e-12: reasons.append('source_activation_bounds')
                if reasons:
                    raise RuntimeError(','.join(reasons))
            guard = ''
            try:
                runner.rollout(u,steps,observer=observe)
            except RuntimeError as error:
                guard = str(error)
            normal,slip = measurements()
            components = dict(terminal_load=float(np.mean(((normal-target)/(.01*bw))**2)),
                              path_load=float(np.mean(path_cost)) if path_cost else 1e6,
                              foot_position=float(np.mean(((d.xpos[feet]-foot0)/.0002)**2)),
                              loaded_slip=float(np.mean((maximum_slip/.005)**2)),
                              root_translation=float(np.mean(((d.qpos[:3]-q0[:3])/.005)**2)),
                              root_rotation=float(np.mean(((d.qpos[3:6]-q0[3:6])/np.deg2rad(2.))**2)))
            if rolling_stance:
                components['source_foot_posture'] = float(np.mean(((d.qpos[foot_posture_address]-q0[foot_posture_address])/foot_posture_limit)**2))
            stance_cost = 'source_foot_posture' if rolling_stance else 'foot_position'
            cost = components['terminal_load']+.1*components['path_load']+sum(components[k] for k in (stance_cost,'loaded_slip','root_translation','root_rotation'))
            if guard: cost += 1e8
            metrics = dict(guard=guard,final_time_s=float(d.time),normal_load_BW=(normal/bw).tolist(),
                           target_normal_load_BW=(target/bw).tolist(),normal_error_BW=((normal-target)/bw).tolist(),
                           maximum_loaded_slip_mps=maximum_slip.tolist(),minimum_visual_mesh_gap_m=minimum_gap.tolist(),
                           maxima=maxima,maximum_knee_extension_deg=np.rad2deg(runner.backend.knee_max).tolist(),components=components)
            if rolling_stance: metrics['contact_consistent_stance'] = contact_maxima
            if store:
                history.append(dict(x=x.tolist(),cost=cost,**metrics))
                if cost<best[0]: best[:] = [cost,x.copy(),metrics]
            if record:
                current_trace[:] = sample
            return cost
        frozen_first = fraction==.01 and unloading==0
        xseed = np.array(accepted_seed['best_group_correction']) if fraction==.01 else previous.copy()
        frozen_extension = extension is not None and any(row['fraction_BW']==fraction for row in extension['stages'] if row['passed'])
        if frozen_extension:
            xseed = np.array(next(row['group_amplitudes'] for row in extension['stages'] if row['fraction_BW']==fraction))
        if replay is not None:
            xseed = np.array(next(row['group_amplitudes'] for row in replay['stages'] if row['fraction_BW']==fraction))
        xseed = np.maximum(low,np.minimum(high,xseed))
        evaluate(xseed)
        if not frozen_first and not frozen_extension and replay is None:
            for j in range(len(groups)):
                for sign in (-1.,1.):
                    probe = xseed.copy(); probe[j] = low[j] if sign<0 else high[j]
                    evaluate(probe)
            constraints = [dict(type='ineq',fun=lambda x:x-low),dict(type='ineq',fun=lambda x:high-x)]
            minimize(evaluate,best[1],method='COBYLA',constraints=constraints,
                     options=dict(rhobeg=.005,tol=1e-5,maxiter=100,catol=1e-10))
        x = best[1]; evaluate(x,store=False,record=True); metrics = best[2]
        complete = abs(float(d.time)-(len(stages)+1)*horizon)<1e-9
        passed = not metrics['guard'] and complete and max(abs(np.array(metrics['normal_error_BW'])))<=.002
        seed_error = None
        if frozen_first and abs(dt-.002)<1e-12:
            prior = np.load(seed_folder/'bounded_states.npz'); index = next(k for k,r in enumerate(seed['rows']) if r['task_feasible'])
            seed_error = {key:float(np.max(abs(value-prior[str(index)+'_'+key]))) for key,value in dict(qpos=d.qpos,qvel=d.qvel,activation=d.act,fiber_state=runner.backend.bank.state).items()}
            if max(seed_error.values())>1e-14: raise RuntimeError('Preserved1% physical replay changed')
        extension_error = None
        if frozen_extension and abs(dt-.002)<1e-12:
            k = int(np.argmin(abs(extension_trace['t']-float(d.time))))
            extension_error = {key:float(np.max(abs(value-extension_trace[key][k]))) for key,value in dict(qpos=d.qpos,qvel=d.qvel,activation=d.act,fiber_state=runner.backend.bank.state,excitation=d.ctrl,source_force_N=runner.backend.bank.force).items()}
            if max(extension_error.values())>1e-14:
                raise RuntimeError('Accepted5% continuation prefix changed')
        trace.extend(current_trace)
        boundary = np.r_[d.qpos,d.qvel,d.act,runner.backend.bank.state,d.ctrl]
        if passed:
            runner.commit_snapshot(); runner.reset()
        boundary_error = float(np.max(abs(boundary-np.r_[d.qpos,d.qvel,d.act,runner.backend.bank.state,d.ctrl])))
        if boundary_error>1e-14: raise RuntimeError('Continuation reset a physical or latent boundary state')
        stages.append(dict(fraction_BW=fraction,passed=bool(passed),frozen_preserved_first_command=frozen_first,
                           frozen_extension_prefix=frozen_extension,extension_prefix_state_error=extension_error,
                           preserved_first_state_error=seed_error,group_amplitudes=x.tolist(),local_lower=low.tolist(),local_upper=high.tolist(),
                           active_local_lower=np.flatnonzero(abs(x-low)<1e-5).tolist(),active_local_upper=np.flatnonzero(abs(x-high)<1e-5).tolist(),
                           boundary_state_error=boundary_error,evaluations=len(history),final=metrics,search_history=history))
        if not passed:
            break
        previous = x; previous_fraction = fraction
    out.mkdir(parents=True,exist_ok=True)
    report = dict(start=str(start_folder),endpoint_phase_reference=str(end_folder),seed=str(seed_folder),unloading_side=('r','l')[unloading],
                  support_fraction=support,timestep_s=dt,horizon_s=horizon,groups=groups,global_group_upper=global_limits.tolist(),
                  source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  rolling_stance=rolling_stance,task_variant='contact_consistent_rolling_stance_v1' if rolling_stance else 'absolute_body_position_hold_v1',
                  extended_from=str(extend_transfer_from) if extend_transfer_from is not None else None,
                  requested_stage_fractions_BW=list(fractions),
                  source_policy='Original source92 and unchanged qualified DS state. Same11family definitions and scalar-amplitude policy both directions; phase-state endpoint differences are not learned L/R gains.',
                  task_policy='Absolute declared-stage BW load targets relative to initial DS; constant segment commands with native50ms filter. C2 load reference enters path cost at0.1weight; terminal force tolerance.002BW. Same absolute physics gates every native step.',
                  stages=stages,native_physics_steps=runner.physics_steps,wall_s=perf_counter()-started,
                  scope='Continuous augmented-state mechanical reference only. Stop at first failed stage; no live SNS,HuntC,source-force or geometry change.')
    if rolling_stance:
        report['rolling_stance_policy'] = dict(tangential_contact_slip_max_mps=.005,visual_mesh_gap_min_m=-.002,
                                              unilateral_normal_tolerance_N=1e-8,friction_cone_tolerance_N=1e-6,
                                              foot_posture_coordinates=foot_posture_names,absolute_foot_posture_max_deg=10.,
                                              posture_bound_provenance='Existing source_force_comparison.full_native_equilibrium non-sagittal/subtalar/MTP teacher bounds +/-10deg. No anatomical range is enlarged; physical model joint limits/geometry remain unchanged.',
                                              position_reference='Actual loaded-contact material-point kinematics. Calcn-body displacement retained as a diagnostic only; removed from hard constraint and soft stance cost in this explicit new variant.',
                                              contact_law='Native2.3.7 pyramidal/elliptic cone,all force/torque components per nativecondim. No added contact force or modified friction.')
    if trace:
        normal_history = np.asarray([row['normal_load_BW'] for row in trace])
        gap_history = np.asarray([row['mesh_gap_m'] for row in trace])
        def longest_duration(mask):
            longest = current = 0
            for value in mask:
                current = current+1 if value else 0
                longest = max(longest,current)
            return longest*dt
        report['actual_unloading_score'] = {}
        for side,label in enumerate(('r','l')):
            force_unloaded = normal_history[:,side]<.02
            mesh_clear = force_unloaded & (gap_history[:,side]>.001)
            force_duration = longest_duration(force_unloaded)
            clear_duration = longest_duration(mesh_clear)
            peak = float(np.max(gap_history[force_unloaded,side])) if np.any(force_unloaded) else None
            report['actual_unloading_score'][label] = dict(minimum_normal_BW=float(np.min(normal_history[:,side])),
                    final_normal_BW=float(normal_history[-1,side]),force_threshold_BW=.02,sustained_duration_threshold_s=.08,
                    longest_force_unloaded_s=force_duration,longest_force_unloaded_mesh_clear_1mm_s=clear_duration,
                    maximum_whole_foot_mesh_gap_m=float(np.max(gap_history[:,side])),
                    peak_whole_foot_gap_while_force_unloaded_m=peak,
                    force_unloading_sustained=bool(force_duration>=.08-1e-10),
                    qualified_clearance_event=bool(clear_duration>=.08-1e-10 and peak is not None and peak>=.01),
                    scope='Separate actual force and whole-foot mesh measurements. Low force alone is not swing. This score does not replace the existing task or physical gates.')
    if replay_folder is not None:
        prior = np.load(replay_folder/'final_continuous_trace.npz')
        now = {key:np.asarray([row[key] for row in trace]) for key in trace[0]}
        keep = np.flatnonzero(prior['t']<=now['t'][-1]+1e-10)
        match = np.searchsorted(now['t'],prior['t'][keep]-1e-10)
        if np.any(match>=len(now['t'])) or np.max(abs(now['t'][match]-prior['t'][keep]))>1e-8:
            raise ValueError('Fresh continuation times do not align')
        report['replay_source'] = str(replay_folder)
        report['replay_maximum_errors'] = {key:float(np.max(abs(now[key][match]-prior[key][keep]))) for key in ('qpos','qvel','qacc','activation','fiber_state','excitation','command','source_force_N','normal_load_BW')}
        report['replay_error_locations'] = {}
        for key in ('qpos','qvel','qacc','activation','fiber_state','source_force_N','normal_load_BW'):
            difference = now[key][match]-prior[key][keep]
            index,column = np.unravel_index(np.argmax(abs(difference)),difference.shape)
            if key in ('qpos','qvel','qacc'):
                coordinate = m.joint(int(m.dof_jntid[column])).name
            elif key=='normal_load_BW':
                coordinate = ('r','l')[column]
            else:
                coordinate = m.actuator(int(column)).name
            report['replay_error_locations'][key] = dict(time_s=float(prior['t'][keep[index]]),coordinate=coordinate,
                                                        reference_value=float(prior[key][keep[index],column]),
                                                        replay_value=float(now[key][match[index],column]),
                                                        signed_error=float(difference[index,column]))
        report['scope'] = 'Fresh saved-command continuous accepted-prefix replay; no refitting. All physical/latent states carry between segments.'
    (out/output_name).write_text(json.dumps(report,indent=2))
    if trace and replay_folder is None:
        np.savez_compressed(out/'final_continuous_trace.npz',group_basis=B,muscles=[m.actuator(j).name for j in range(m.nu)],parameters=runner.backend.parameters,
                            **{key:np.asarray([row[key] for row in trace]) for key in trace[0]})
    print(json.dumps({k:v for k,v in report.items() if k!='stages'},indent=2))
    print(json.dumps([{k:v for k,v in row.items() if k!='search_history'} for row in stages],indent=2))


def audit_transfer_geometry(model_folder,start_folder,folder,out):
    """Read-only FK decomposition of a stopped continuous-transfer reference."""
    from scipy.spatial.transform import Rotation
    m = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
    d = mujoco.MjData(m); start = np.load(start_folder/'source_recruitment.npz')
    z = np.load(folder/'final_continuous_trace.npz')
    d.qpos[:] = start['qpos']; mujoco.mj_fwdPosition(m,d)
    feet = [m.body('calcn_'+side).id for side in ('r','l')]
    p0 = d.xpos[feet].copy(); r0 = d.xmat[feet].reshape(2,3,3).copy()
    maximum = dict(displacement_m=-1.); final = []
    for k,q in enumerate(z['qpos']):
        d.qpos[:] = q; mujoco.mj_fwdPosition(m,d)
        difference = d.xpos[feet]-p0
        side,axis = np.unravel_index(np.argmax(abs(difference)),difference.shape)
        if abs(difference[side,axis])>maximum['displacement_m']:
            maximum = dict(time_s=float(z['t'][k]),side=('r','l')[side],axis=('anterior_x','left_y','superior_z')[axis],
                           displacement_m=float(abs(difference[side,axis])),signed_displacement_xyz_m=difference[side].tolist(),
                           world_rotation_vector_deg=np.rad2deg(Rotation.from_matrix(d.xmat[feet[side]].reshape(3,3)@r0[side].T).as_rotvec()).tolist())
        if k==len(z['qpos'])-1:
            final = [dict(side=('r','l')[s],displacement_xyz_m=difference[s].tolist(),world_rotation_vector_deg=np.rad2deg(Rotation.from_matrix(d.xmat[b].reshape(3,3)@r0[s].T).as_rotvec()).tolist()) for s,b in enumerate(feet)]
    result = dict(source=str(folder),frame='Captured human native world:xanterior,yleft,zsuperior; no robot-frame reinterpretation.',
                  maximum_absolute_calcn_body_position_error=maximum,final_foot_body_changes=final,
                  interpretation='Body-origin displacement includes foot rotation/rocking and compliant deformation; it is not identical to measured loaded-contact slip. This audit does not relax either constraint.')
    out.mkdir(parents=True,exist_ok=True); target = out/'TRANSFER_GEOMETRY_AUDIT.json'
    if target.exists(): raise ValueError('Preserve geometry audit')
    target.write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))


def audit_transfer_authority(model_folder,start_folder,folder,out):
    """Read-only limits, delivered recruitment, and endpoint contact ownership."""
    report = json.loads((folder/'CONTINUOUS_GROUP_TRANSFER.json').read_text())
    start = np.load(start_folder/'source_recruitment.npz'); z = np.load(folder/'final_continuous_trace.npz')
    m = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb')); d = mujoco.MjData(m)
    d.qpos[:] = z['qpos'][-1]; d.qvel[:] = z['qvel'][-1]; d.act[:] = z['activation'][-1]
    d.qacc_warmstart[:] = z['qacc'][-1]
    d.xfrc_applied[m.body('pelvis').id,2] = report['support_fraction']*float(m.body_mass.sum()*9.81)
    m.actuator_gaintype[:] = mujoco.mjtGain.mjGAIN_FIXED; m.actuator_gainprm[:] = 0.
    m.actuator_biastype[:] = mujoco.mjtBias.mjBIAS_AFFINE; m.actuator_biasprm[:] = 0.
    m.actuator_biasprm[:,0] = -z['source_force_N'][-1]; mujoco.mj_forward(m,d)
    probe = NativeContactKinematics(m,d); measured = probe.measure(); bw = float(m.body_mass.sum()*9.81)
    normal_by_body = {name:0. for side in ('r','l') for name in ('calcn_'+side,'toes_'+side)}
    for k in range(d.ncon):
        c = d.contact[k]
        if c.efc_address<0: continue
        mujoco.mj_contactForce(m,d,k,probe.force)
        for body in (int(m.geom_bodyid[c.geom1]),int(m.geom_bodyid[c.geom2])):
            name = m.body(body).name
            if name in normal_by_body: normal_by_body[name] += float(probe.force[0])
    last = report['stages'][-1]; origin = start['activation']; names = z['muscles'].tolist(); limits = []
    for j in last['active_local_upper']:
        basis = z['group_basis'][:,j]; positive = basis>0.; negative = basis<0.
        ceiling = min(float(np.min((.98-origin[positive])/basis[positive])) if np.any(positive) else np.inf,
                      float(np.min((origin[negative]-.01)/(-basis[negative]))) if np.any(negative) else np.inf)
        limits.append(dict(group=report['groups'][j],amplitude=last['group_amplitudes'][j],
                           local_upper=last['local_upper'][j],endpoint_derived_global_upper=report['global_group_upper'][j],
                           physical_command_upper_along_same_direction=ceiling,
                           muscles=[dict(name=n,initial_command=float(origin[i]),basis_coefficient=float(basis[i]),
                                         final_command=float(z['command'][-1,i]),delivered_activation=float(z['activation'][-1,i]),
                                         delivered_force_N=float(z['source_force_N'][-1,i])) for i,n in enumerate(names) if basis[i]!=0.]))
    result = dict(source=str(folder),last_target_BW=last['fraction_BW'],last_target_passed=last['passed'],last_guard=last['final']['guard'],
                  actual_native_root_delta={m.joint(int(m.dof_jntid[k])).name:float(d.qpos[k]-start['qpos'][k]) for k in range(6)},
                  baseline_anatomical_angles_deg={name:float(np.rad2deg(start['qpos'][m.joint(name).qposadr[0]])) for name in JOINTS},
                  active_recruitment_bounds=limits,endpoint_contact_normal_N_by_body=normal_by_body,
                  reconstructed_normal_error_BW=float(np.max(abs(measured['normal_N']/bw-z['normal_load_BW'][-1]))),
                  scope='Scratch instantaneous force/contact reconstruction checked against saved normal loads. No forward trial or parameter fit. Recruitment-direction/local-box limits are not global muscle-strength limits.')
    out.mkdir(parents=True,exist_ok=True); target = out/'TRANSFER_AUTHORITY_AUDIT.json'
    if target.exists(): raise ValueError('Preserve authority audit')
    target.write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))


def prepare_stationary_transfer(model_folder,start_folder,folder,out,fraction=.125):
    """Export an achieved loaded pose as an explicit static inverse seed."""
    report = json.loads((folder/'CONTINUOUS_GROUP_TRANSFER.json').read_text())
    stage = next(row for row in report['stages'] if row['fraction_BW']==fraction and row['passed'])
    z = np.load(folder/'final_continuous_trace.npz'); origin = np.load(start_folder/'source_recruitment.npz')
    index = int(np.argmin(abs(z['t']-stage['final']['final_time_s'])))
    if abs(z['t'][index]-stage['final']['final_time_s'])>1e-8: raise ValueError('Accepted static seed time missing')
    m = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb')); d = mujoco.MjData(m)
    d.qpos[:] = z['qpos'][index]; mujoco.mj_fwdPosition(m,d)
    parameters = z['parameters']; alpha = z['activation'][index]
    bank = DynamicThelenBank(parameters,float(m.opt.timestep)); bank.initialize(d.actuator_length,alpha)
    achieved = z['normal_load_BW'][index]
    total = float(m.body_subtreemass[m.body('pelvis').id]/m.body_mass.sum()-report['support_fraction'])
    target = achieved/achieved.sum()*total
    reference = dict(source_transfer=str(folder),fraction_BW=fraction,source_time_s=float(z['t'][index]),
                     achieved_dynamic_normal_BW=achieved.tolist(),normal_target_BW=target.tolist(),
                     root_anchor_native=origin['qpos'][:6].tolist(),
                     normalization='Preserve achieved R/L split ratio; static total normal equals floating-body weight minus unchanged prescribed support. Vertical rig stiffness is zero; no force is added.',
                     reference_policy='Actual loaded rolling foot pose and105 native coordinates; no follower projection in seed. New static activation states are not learned neural gains.')
    manifest = dict(format='aarl_rolling_load_stationary_seed_v1',support_fraction=report['support_fraction'],stationary_load_split_reference=reference)
    out.mkdir(parents=True,exist_ok=True)
    if (out/'source_recruitment.npz').exists(): raise ValueError('Preserve stationary seed')
    np.savez_compressed(out/'source_recruitment.npz',qpos=d.qpos.copy(),activation=alpha,source_force=bank.force,
                        parameters=parameters,muscles=z['muscles'],force_channel_enabled=np.ones(m.nu,bool))
    (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    (out/'SOURCE_NATIVE_STANCE.json').write_text(json.dumps(dict(source_pose=str(out),accepted_static_source_recruitment=False,accepted_stationary_endpoint=False,scope='Prepared seed only; static inverse and five-second validation pending.'),indent=2))
    print(json.dumps(manifest,indent=2))


def prepare_manifold_seed(model_folder,start_folder,end_folder,out,fraction=.01):
    """Small load-split initial guess between frozen equilibria,not acceptance."""
    from source_force_comparison import load_recruitment
    if not 0.<fraction<=.05:raise ValueError('Declared small-manifold task is limited to5% BW')
    if (out/'source_recruitment.npz').exists():raise ValueError('Preserve existing manifold seed')
    source=[json.loads((f/'SOURCE_NATIVE_STANCE.json').read_text()) for f in (start_folder,end_folder)]
    if not all(row.get('accepted_stationary_endpoint',False) for row in source):raise ValueError('Accepted stationary endpoints required')
    a,b=load_recruitment(start_folder),load_recruitment(end_folder)
    if not np.array_equal(a['parameters'],b['parameters']) or not np.array_equal(a['muscles'],b['muscles']):raise ValueError('Endpoint source parameters differ')
    endcfg=json.loads((end_folder/'MANIFEST.json').read_text());support=endcfg['support_fraction']
    endnormal=np.array(endcfg['stationary_load_split_reference']['normal_target_BW'])
    normal=np.load(start_folder/'stationarity_sparse.npz')['normal_side_BW'][-50:].mean(axis=0)
    m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m)
    total=float(m.body_subtreemass[m.body('pelvis').id]/m.body_mass.sum()-support);normal*=total/normal.sum()
    direction=np.sign(endnormal[0]-normal[0]);ratio=fraction/abs(endnormal[0]-normal[0])
    q=a['qpos']+ratio*(b['qpos']-a['qpos']);alpha=a['activation']+ratio*(b['activation']-a['activation'])
    d.qpos[:]=q;mujoco.mj_fwdPosition(m,d);bank=DynamicThelenBank(a['parameters'],float(m.opt.timestep));bank.initialize(d.actuator_length,alpha)
    reference=dict(source_transfer=str(end_folder),fraction_BW=fraction,normal_target_BW=(normal+direction*fraction*np.array([1.,-1.])).tolist(),
                   root_anchor_native=a['qpos'][:6].tolist(),normalization='Exact declared BW transfer from validated DS,stationary total normal retained; no additional force.',
                   reference_policy='Linear full-native pose/activation initialization only; all105 equilibrium equations and5-second stationarity remain required.')
    manifest=dict(format='aarl_small_stationary_manifold_seed_v1',support_fraction=support,stationary_load_split_reference=reference,
                  frozen_start=str(start_folder),frozen_end=str(end_folder),initialization_ratio=float(ratio),source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out.mkdir(parents=True,exist_ok=True);np.savez_compressed(out/'source_recruitment.npz',qpos=q,activation=alpha,source_force=bank.force,parameters=a['parameters'],muscles=a['muscles'],force_channel_enabled=a['force_channel_enabled'])
    (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));(out/'SOURCE_NATIVE_STANCE.json').write_text(json.dumps(dict(source_pose=str(out),accepted_static_source_recruitment=False,accepted_stationary_endpoint=False),indent=2))
    print(json.dumps(manifest,indent=2))


def stationary_source_capacity(model_folder,folder,out):
    """Necessary source-force recruitment at fixed loaded geometry,not equilibrium."""
    from scipy.optimize import lsq_linear
    from contact_standing_solver import _external
    saved=np.load(folder/'source_recruitment.npz');cfg=json.loads((folder/'MANIFEST.json').read_text())
    m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));d=mujoco.MjData(m)
    d.qpos[:]=saved['qpos'];d.qvel[:]=0.;d.qacc[:]=0.;d.xfrc_applied[m.body('pelvis').id,2]=cfg['support_fraction']*float(m.body_mass.sum()*9.81)
    mujoco.mj_inverse(m,d);ar=CoupledArms(m);arms=ar.update(d).copy();required=ar.E.T@(d.qfrc_inverse-_external(m,d))
    bank=DynamicThelenBank(saved['parameters'],float(m.opt.timestep));alpha=saved['activation'].copy();initial=alpha.copy();history=[]
    def force(a):
        bank.initialize(d.actuator_length,a)
        if np.any(bank.bounds):raise RuntimeError('Static source fiber bound during capacity audit')
        return bank.force.copy()
    for iteration in range(8):
        current=force(alpha);lo=np.maximum(.01,alpha-.005);hi=np.minimum(.98,alpha+.005)
        gain=(force(hi)-force(lo))/(hi-lo)
        if min(gain)<-1e-5:raise RuntimeError('Nonmonotone source static activation response')
        A=arms*np.maximum(gain,0.)[None,:]
        fit=lsq_linear(np.vstack([A,.0001*np.eye(m.nu)]),np.r_[required-arms@current+A@alpha,.0001*initial],bounds=(np.full(m.nu,.01),np.full(m.nu,.98)),method='bvls',tol=1e-10,max_iter=200)
        alpha=fit.x;residual=arms@force(alpha)-required
        history.append(dict(iteration=iteration,anatomical_torque_residual_max_Nm=float(np.max(abs(residual))),rms_Nm=float(np.sqrt(np.mean(residual**2)))))
        if max(abs(residual))<.001:break
    F=force(alpha);raw=-d.actuator_moment.T@F-(d.qfrc_inverse-_external(m,d))
    result=dict(source=str(folder),anatomical_joint_names=list(JOINTS),anatomical_required_torque_Nm=required.tolist(),
                anatomical_residual_Nm=(arms@F-required).tolist(),full105_raw_residual_max_N_or_Nm=float(np.max(abs(raw))),
                follower_force_residual_max_N=float(np.max(abs(raw[ar.dependent]))),history=history,
                active_command_floor=[m.actuator(j).name for j in np.flatnonzero(alpha<.01001)],active_command_ceiling=[m.actuator(j).name for j in np.flatnonzero(alpha>.97999)],
                scope='Necessary equality-aware anatomical torque capacity only. All105 raw/follower residuals are reported and remain required; no projected-cancellation acceptance. Activations are source states,not neural gains; all strength parameters remain tied.')
    out.mkdir(parents=True,exist_ok=True)
    if (out/'SOURCE_CAPACITY.json').exists():raise ValueError('Preserve source capacity audit')
    (out/'SOURCE_CAPACITY.json').write_text(json.dumps(result,indent=2))
    data={k:saved[k] for k in saved.files};data['activation']=alpha;data['source_force']=F
    if 'fit_variables' in data:
        x=data['fit_variables'].copy();x[m.nv:m.nv+m.nu]=alpha;data['fit_variables']=x
    np.savez_compressed(out/'source_recruitment.npz',**data)
    cfg['source_capacity_reference']=str(folder);cfg['source_capacity_scope']=result['scope'];(out/'MANIFEST.json').write_text(json.dumps(cfg,indent=2))
    (out/'SOURCE_NATIVE_STANCE.json').write_text(json.dumps(dict(source_pose=str(out),accepted_static_source_recruitment=False,accepted_stationary_endpoint=False,scope=result['scope']),indent=2))
    print(json.dumps(result,indent=2))


def shoot_stationary_manifold(model_folder,start_folder,target_folders,spec_path,out,budget=24,physics_dt=None,terminal_window=False,replay_candidate=None,full_family_audit=False,fd_step=.0001,derivative_free_task=False,proposal_audit=None):
    """Bounded paired nonlinear shooting with one shared source-family waveform."""
    from time import perf_counter
    from scipy.optimize import least_squares
    from contact_standing_solver import _external
    from source_force_comparison import load_recruitment
    if budget<(3 if replay_candidate is not None else 10):raise ValueError('Budget must include initialization,audit,and reserved final replay')
    if not np.isfinite(fd_step) or not 0.<fd_step<=.05:raise ValueError('Finite difference step must lie in(0,.05]')
    if derivative_free_task and (budget>24 or physics_dt!=.001 or not terminal_window or full_family_audit or replay_candidate is not None):raise ValueError('Task search requires <=24 paired evaluations,1ms,terminal window,and no gradient/candidate mode')
    if (out/'MANIFOLD_SHOOTING.json').exists():raise ValueError('Preserve completed shooting')
    out.mkdir(parents=True,exist_ok=True);started=perf_counter();spec=json.loads(spec_path.read_text())
    start=load_recruitment(start_folder);targets=[load_recruitment(f) for f in target_folders]
    for folder in [start_folder]+list(target_folders):
        if not json.loads((folder/'SOURCE_NATIVE_STANCE.json').read_text()).get('accepted_stationary_endpoint',False):raise ValueError('Validated stationary references required')
    if any(not np.array_equal(t['parameters'],start['parameters']) or not np.array_equal(t['muscles'],start['muscles']) for t in targets):raise ValueError('Source endpoint parameters differ')
    names=start['muscles'].tolist();families=sorted({n[:-2] for n in names});p=start['parameters'];a0=start['activation']
    if not start['force_channel_enabled'].all() or any(not t['force_channel_enabled'].all() for t in targets):raise ValueError('Explicit original92-source policy required')
    address=np.array([v['qpos_address'] for v in spec['joints'].values()]);low=np.array([v['native_bounds_rad'][0] for v in spec['joints'].values()]);high=np.array([v['native_bounds_rad'][1] for v in spec['joints'].values()])
    reference=np.load(spec['source_target']);vmax=np.max(abs(reference['qvel'][:,address]),axis=0);amax=np.max(abs(reference['qacc'][:,address]),axis=0)
    warm,ramp,hold,filter_tau=.5,1.812,5.,.05
    branches=[];torque_maps=[]
    for folder,target in zip(target_folders,targets):
        m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
        if m.nu!=92 or [m.actuator(i).name for i in range(m.nu)]!=names:raise ValueError('Captured original92 actuator order required')
        if physics_dt is not None:m.opt.timestep=physics_dt
        r=AugmentedTransition(m,p,filter_tau);m,d,b=r.model,r.data,r.backend;dt=float(m.opt.timestep)
        d.qpos[:]=start['qpos'];d.qvel[:]=0.;d.act[:]=a0;d.ctrl[:]=a0
        bw=float(m.body_mass.sum()*9.81);d.xfrc_applied[m.body('pelvis').id,2]=.6*bw;mujoco.mj_forward(m,d)
        b.bank.initialize(d.actuator_length,a0);m.actuator_biasprm[:,0]=-b.bank.force;mujoco.mj_forward(m,d);r.commit_snapshot()
        ar=CoupledArms(m);probe=mujoco.MjData(m);probe.qpos[:]=target['qpos'];probe.act[:]=target['activation'];mujoco.mj_forward(m,probe);ar.update(probe)
        bank=DynamicThelenBank(p,dt);al=np.maximum(.01,target['activation']-.0001);ah=np.minimum(.98,target['activation']+.0001)
        bank.initialize(probe.actuator_length,al);fl=bank.force.copy();bank.initialize(probe.actuator_length,ah);gain=(bank.force-fl)/(ah-al)
        B=np.array([[target['activation'][k]-a0[k] if n[:-2]==family else 0. for k,n in enumerate(names)] for family in families])
        torque_maps.append((ar.E.T@(-probe.actuator_moment.T*gain[None,:]))@B.T)
        bank.initialize(probe.actuator_length,target['activation']);target_fiber=bank.state.copy();geometry=[]
        for g in range(m.ngeom):
            name=m.body(int(m.geom_bodyid[g])).name
            if m.geom_type[g]==mujoco.mjtGeom.mjGEOM_MESH and name in ('calcn_r','toes_r','calcn_l','toes_l'):
                mesh=int(m.geom_dataid[g]);first=int(m.mesh_vertadr[mesh]);number=int(m.mesh_vertnum[mesh]);geometry.append((0 if name.endswith('_r') else 1,g,m.mesh_vert[first:first+number].copy()))
        cfg=json.loads((folder/'MANIFEST.json').read_text());normal=np.array(cfg['stationary_load_split_reference']['normal_target_BW'])
        branches.append(dict(runner=r,arms=ar,probe=probe,static_bank=bank,target=target,target_fiber=target_fiber,basis=B,
                             normal_target=normal,bw=bw,geometry=geometry,contact=NativeContactKinematics(m,d),folder=folder))
    task_map=np.vstack(torque_maps);row_norm=np.linalg.norm(task_map,axis=1);task_normalized=task_map/np.maximum(1e-12,row_norm[:,None])
    task_labels=None
    if derivative_free_task:
        task_labels=['hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','subtalar_angle','mtp_angle','lumbar_extension','lumbar_bending','lumbar_rotation']
        wanted=np.zeros((len(branches)*len(JOINTS),len(task_labels)))
        for phase in range(2):
            for j in range(7):wanted[phase*17+j,j]=1.-2*phase;wanted[phase*17+7+j,j]=-1.+2*phase
            for j in range(3):wanted[phase*17+14+j,7+j]=1. if j==0 else 1.-2*phase
        # Explicit normalized anatomical task directions,not global torque modes.
        modes=np.linalg.solve(task_normalized.T@task_normalized+1e-6*np.eye(len(families)),task_normalized.T@wanted)
        norms=np.linalg.norm(modes,axis=0)
        if np.any(norms<1e-12):raise ValueError('An anatomical task direction has no source-family authority')
        modes/=norms;singular=np.linalg.svd(task_normalized,compute_uv=False)
    else:
        _,singular,Vh=np.linalg.svd(task_normalized if full_family_audit else task_map,full_matrices=full_family_audit)
        modes=Vh.T if full_family_audit else Vh[:6].T
    nmode=modes.shape[1]
    if full_family_audit and budget<nmode+7:raise ValueError('Full family audit needs baseline,reset,all columns,scale checks,and final replay')
    manifest=dict(format='aarl_shared_stationary_manifold_shooting_v1',source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  start=str(start_folder),targets=[str(f) for f in target_folders],support_fraction=.6,warm_s=warm,ramp_s=ramp,hold_s=hold,filter_tau_s=filter_tau,
                  budget_paired_rollouts=budget,families=families,source_torque_singular_values=singular.tolist(),shared_modes=nmode,
                  full_family_audit=full_family_audit,finite_difference_step=fd_step,task_row_norms=row_norm.tolist(),
                  retained_squared_row_coverage=np.divide(np.linalg.norm(task_map@modes,axis=1)**2,np.maximum(1e-30,row_norm**2)).tolist(),
                  rule='Common coefficients across homologous muscles AND both phases. Full task-row-normalized orthogonal family basis including nullspace; no truncation.' if full_family_audit else 'Common coefficients across homologous muscles AND both phases; six global torque directions only.',
                  waveform='alpha0+s*delta+16*s^2*(1-s)^2*sum(c_family*delta_family); s is quintic C2 phase. Native activation,fiber,filter,and contacts propagate without trajectory resets.',
                  terminal_window_objective=terminal_window,candidate_file=str(replay_candidate) if replay_candidate is not None else None,
                  candidate_sha256=hashlib.sha256(replay_candidate.read_bytes()).hexdigest() if replay_candidate is not None else None,
                  scope='Small loaded mechanical reference only; no SNS fit,Hunt C promotion,or swing claim.')
    if derivative_free_task:
        projection=modes@np.linalg.pinv(modes)
        manifest.update(derivative_free_task=True,anatomical_task_labels=task_labels,optimization_dimension=2*nmode,
            rule='Ten explicit normalized anatomical directions retaining ankle,subtalar,and MTP actions. The same family coefficients apply to both limbs and both phases; no global torque SVD truncation.',
            waveform='Existing C2 ramp plus pulse*family_c0+pulse*(2*s-1)*family_c1; both corrections vanish at immutable endpoints.',
            trust_region='Each mode coefficient <=0.1; maximum per-family sum of absolute temporal coefficients <=0.1.',
            anatomical_row_projection_coverage=(np.linalg.norm(task_normalized@projection,axis=1)**2).tolist(),
            retained_squared_row_coverage=(np.linalg.norm(task_normalized@projection,axis=1)**2).tolist(),
            proposal_audit=str(proposal_audit) if proposal_audit else None,
            score='Maximum and mean log1p squared actual acceptance-threshold ratios across both phases; physical guards carry a separate dominant penalty. No force/velocity gradient fit.')
    (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    if derivative_free_task:
        distal=[phase*17+leg*7+j for phase in range(2) for leg in range(2) for j in (4,5,6)]
        coverage=np.array(manifest['anatomical_row_projection_coverage'])
        manifest['anatomical_basis_gate']=dict(required_distal_squared_row_coverage=.95,measured_minimum=float(min(coverage[distal])),accepted=bool(min(coverage[distal])>=.95))
        (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
        if not manifest['anatomical_basis_gate']['accepted']:raise ValueError('Reduced anatomical basis fails distal authority preflight; do not launch shooting')
    history=[];best=[np.inf,None,None];cache={};calls=0
    def branch(branch,coefficients,record=False):
        r=branch['runner'];r.reset();m,d,b=r.model,r.data,r.backend;ar=branch['arms'];dt=float(m.opt.timestep)
        target=branch['target'];n=round((warm+ramp+hold)/dt);late_start=n-round(.5/dt);qgoal=target['qpos'];normal=np.zeros(2)
        if full_family_audit:
            import zlib
            if 'branch_buffers' not in branch:branch['branch_buffers']=(np.empty((n,19)),np.empty((n,len(names)),np.uint8),np.empty((n,len(names)),np.uint8))
            phase_history,muscle_history,wrap_history=branch['branch_buffers'];mtp=np.array([int(m.joint('mtp_angle_'+s).qposadr[0]) for s in ('r','l')]);tendon_ids=m.actuator_trnid[:,0]
        driver=np.r_[np.arange(6),ar.addresses];qscale=np.r_[[.0005]*3,[np.deg2rad(.25)]*(len(driver)-3)];vscale=np.r_[[.0001]*3,[np.deg2rad(.05)]*(len(driver)-3)]
        max_slip=0.;max_cone=0.;max_root=0.;max_rotation=0.;min_gap=np.full(2,np.inf);max_acc=0.;max_speed=0.;guard='';late=[];hold_q=None;trace=[];hold_root_drift=0.;hold_joint_drift=0.;hold_rotation_drift=0.
        checksum=hashlib.sha256() if derivative_free_task else None
        correction=coefficients@branch['basis']
        for i in range(n):
            phase=np.clip((d.time-warm)/ramp,0.,1.);s=phase**3*(10.+phase*(-15.+6.*phase));pulse=16*s*s*(1-s)*(1-s)
            command=a0+s*(target['activation']-a0)+pulse*(correction[0]+(2*s-1)*correction[1] if derivative_free_task else correction)
            if command.min()<.01-1e-12 or command.max()>.98+1e-12:guard='Source command bounds';break
            r.excitation+=r.beta*(command-r.excitation);d.ctrl[:]=r.excitation
            try:b.step();r.physics_steps+=1
            except RuntimeError as error:guard=str(error);break
            if checksum is not None:checksum.update(np.r_[d.qpos,d.qvel,d.act,b.bank.state,b.bank.force,d.ctrl].tobytes())
            contacts=branch['contact'].measure(full_family_audit);normal=contacts['normal_N']/branch['bw'];max_slip=max(max_slip,float(max(contacts['slip_mps'])));max_cone=max(max_cone,contacts['maximum_friction_cone_violation_N'])
            if full_family_audit:
                length=np.hypot(b.bank.state,p[:,3]);velocity=(b.bank.state-b.old_state)*b.bank.state/(length*dt*p[:,1]);c=.25+.75*np.maximum(.01,b.activation);fvz=velocity/p[:,8];k=(2.+2./p[:,6])/(p[:,7]-1.);edge=.95*p[:,7];threshold=c*(edge-1.)/(k*(p[:,7]-edge))
                fv=np.select([fvz<=-c,velocity<=0.,fvz<=threshold],[0,1,2],default=3).astype(np.uint8);strain=(d.actuator_length-b.bank.state)/p[:,2]-1.;toe=99.*np.exp(3.)/(166.*np.exp(3.)-67.)*p[:,4]
                tendon=np.select([strain<=0.,strain<toe],[0,1],default=2).astype(np.uint8)
                muscle_history[i]=fv+4*(length>p[:,1])+8*tendon+32*(d.ctrl>b.activation)
                wrap_history[i]=d.ten_wrapnum[tendon_ids]
                phase_history[i]=np.r_[d.time,d.qpos[mtp],d.qvel[mtp],d.qacc[mtp],contacts['normal_by_body_N'],contacts['contact_count_by_body'],contacts['signature_u32'],int(d.solver_iter),d.nefc,zlib.crc32(wrap_history[i].tobytes())]
            gap=np.full(2,np.inf)
            for side,g,vertices in branch['geometry']:gap[side]=min(gap[side],float(np.min(vertices@d.geom_xmat[g].reshape(3,3)[2]+d.geom_xpos[g,2])))
            min_gap=np.minimum(min_gap,gap);root=float(max(abs(d.qpos[:3]-start['qpos'][:3])));rotation=float(max(abs(d.qpos[3:6]-start['qpos'][3:6])))
            max_root=max(max_root,root);max_rotation=max(max_rotation,rotation);max_acc=max(max_acc,float(max(abs(d.qacc[address])/amax)));max_speed=max(max_speed,float(max(abs(d.qvel[address])/vmax)))
            if np.any(d.qpos[address]<low-1e-10) or np.any(d.qpos[address]>high+1e-10):guard='Frozen source ROM'
            elif root>.005 or rotation>np.deg2rad(2.):guard='Root envelope'
            elif max_acc>1.+1e-10 or max_speed>1.+1e-10:guard='Source rate envelope'
            elif max_slip>.005:guard='Loaded contact slip'
            elif min_gap.min()<-.002:guard='Visual mesh penetration'
            elif max_cone>1e-6 or contacts['minimum_normal_N']< -1e-8:guard='Native contact cone'
            other=np.array([int(m.joint(int(j)).qposadr[0]) for j in ar.jids if m.joint(int(j)).name not in spec['joints']])
            if np.any(abs(d.qpos[other])>np.deg2rad(10.)+1e-10):guard=guard or 'Existing non-sagittal source bounds'
            if hold_q is None and d.time>=warm+ramp-1e-10:hold_q=d.qpos.copy()
            if hold_q is not None:
                hold_root_drift=max(hold_root_drift,float(max(abs(d.qpos[:3]-hold_q[:3]))));hold_joint_drift=max(hold_joint_drift,float(max(abs(d.qpos[ar.addresses]-hold_q[ar.addresses]))))
                hold_rotation_drift=max(hold_rotation_drift,float(max(abs(d.qpos[3:6]-hold_q[3:6]))))
            if i>=late_start:
                length=np.hypot(b.bank.state,p[:,3]);fv=(b.bank.state-b.old_state)*b.bank.state/(length*dt*p[:,1])
                if terminal_window:
                    probe=branch['probe'];probe.qpos[:]=d.qpos;probe.qvel[:]=0.;probe.qacc[:]=0.;probe.act[:]=d.act;probe.xfrc_applied[:]=d.xfrc_applied;mujoco.mj_inverse(m,probe)
                    required=probe.qfrc_inverse-_external(m,probe);branch['static_bank'].initialize(probe.actuator_length,d.act)
                    window_dynamic=-probe.actuator_moment.T@b.bank.force-required;window_static=-probe.actuator_moment.T@branch['static_bank'].force-required
                else:window_dynamic=window_static=np.zeros(m.nv)
                late.append((d.qpos.copy(),d.qvel.copy(),b.bank.force.copy(),fv.copy(),normal.copy(),b.bank.state.copy(),window_dynamic.copy(),window_static.copy()))
            if record and (i%max(1,round(.01/dt))==0 or i==n-1):trace.append(np.r_[d.time,d.qpos,d.qvel,d.act,b.bank.state,d.ctrl,normal,gap])
            if guard:break
        probe=branch['probe'];probe.qpos[:]=d.qpos;probe.qvel[:]=0.;probe.qacc[:]=0.;probe.act[:]=d.act;probe.xfrc_applied[:]=d.xfrc_applied;mujoco.mj_inverse(m,probe)
        required=probe.qfrc_inverse-_external(m,probe);branch['static_bank'].initialize(probe.actuator_length,d.act)
        raw=-probe.actuator_moment.T@b.bank.force-required;static_raw=-probe.actuator_moment.T@branch['static_bank'].force-required
        residual=np.r_[(d.qpos[driver]-qgoal[driver])/qscale,d.qvel[driver]/vscale,(b.bank.state-branch['target_fiber'])/(.001*p[:,1]),raw/.05,static_raw/.05,(normal-branch['normal_target'])/.002,0. if not guard else 100.*(n-i)/n]
        point_residual=residual.copy()
        if terminal_window and late:
            qlate=np.array([row[0][driver] for row in late]);vlate=np.array([row[1][driver] for row in late]);flate=np.array([row[5] for row in late]);dynamic=np.array([row[6] for row in late]);static=np.array([row[7] for row in late]);nlate=np.array([row[4] for row in late])
            residual=np.r_[(qlate.mean(axis=0)-qgoal[driver])/qscale,vlate.mean(axis=0)/vscale,(flate.mean(axis=0)-branch['target_fiber'])/(.001*p[:,1]),dynamic.mean(axis=0)/.05,static.mean(axis=0)/.05,(nlate.mean(axis=0)-branch['normal_target'])/.002,
                           vlate.std(axis=0)/vscale,dynamic.std(axis=0)/.05,static.std(axis=0)/.05,0. if not guard else 100.*(n-i)/n]
        elif terminal_window:
            residual=np.r_[point_residual[:-1],np.zeros(len(driver)+2*m.nv),point_residual[-1]]
        late_speed=max((float(max(abs(row[1][ar.rows]))) for row in late),default=np.inf);late_root_speed=max((float(max(abs(row[1][:3]))) for row in late),default=np.inf)
        fiber_speed=max((float(max(abs(row[3]))) for row in late),default=np.inf);force_variation=float(np.max(np.ptp(np.array([row[2] for row in late]),axis=0)/p[:,0])) if late else np.inf
        normal_error=max((float(max(abs(row[4]-branch['normal_target']))) for row in late),default=np.inf)
        hold_drift=hold_root_drift if hold_q is not None else np.inf
        pose_error=float(max(abs(d.qpos[ar.addresses]-qgoal[ar.addresses])));root_pose_error=float(max(abs(d.qpos[:3]-qgoal[:3])))
        late_root_angular_speed=max((float(max(abs(row[1][3:6]))) for row in late),default=np.inf);root_rotation_error=float(max(abs(d.qpos[3:6]-qgoal[3:6])))
        endpoint_activation_error=float(max(np.max(abs(d.act-target['activation'])),np.max(abs(d.ctrl-target['activation']))))
        accepted=bool(not guard and i==n-1 and normal_error<=.002 and max(abs(raw))<=.05 and max(abs(static_raw))<=.05 and late_speed<=np.deg2rad(.05) and late_root_speed<=.0001 and late_root_angular_speed<=np.deg2rad(.05) and fiber_speed<=.001 and force_variation<=1e-4 and pose_error<=np.deg2rad(.25) and root_pose_error<=.0005 and root_rotation_error<=np.deg2rad(.25) and hold_drift<=.0005 and hold_joint_drift<=np.deg2rad(.25) and hold_rotation_drift<=np.deg2rad(.25) and endpoint_activation_error<=1e-12)
        metrics=dict(accepted=accepted,time_s=float(d.time),guard=guard,maximum_loaded_slip_mps=max_slip,maximum_root_translation_m=max_root,maximum_root_rotation_deg=float(np.rad2deg(max_rotation)),maximum_friction_cone_violation_N=max_cone,minimum_mesh_gap_m=min_gap.tolist(),source_acceleration_ratio=max_acc,source_velocity_ratio=max_speed,
                     normal_BW=normal.tolist(),normal_target_BW=branch['normal_target'].tolist(),late_normal_error_BW=normal_error,final_dynamic_force_static_residual_max=float(max(abs(raw))),final_static_residual_max=float(max(abs(static_raw))),late_joint_speed_deg_s=float(np.rad2deg(late_speed)),late_root_speed_mps=late_root_speed,late_fiber_speed_l0_s=fiber_speed,late_force_variation_F0=force_variation,final_joint_error_deg=float(np.rad2deg(pose_error)),final_root_error_m=root_pose_error,terminal_hold_root_drift_m=hold_drift,maximum_knee_extension_deg=np.rad2deg(b.knee_max).tolist())
        state=np.r_[d.qpos,d.qvel,d.act,b.bank.state,d.ctrl]
        metrics['terminal_hold_joint_drift_deg']=float(np.rad2deg(hold_joint_drift))
        metrics.update(terminal_hold_root_rotation_drift_deg=float(np.rad2deg(hold_rotation_drift)),late_root_angular_speed_deg_s=float(np.rad2deg(late_root_angular_speed)),final_root_rotation_error_deg=float(np.rad2deg(root_rotation_error)),final_activation_excitation_error=endpoint_activation_error)
        if checksum is not None:metrics['physical_state_force_sha256']=checksum.hexdigest()
        if full_family_audit:
            count=i+1;current=(phase_history[:count],muscle_history[:count],wrap_history[:count]);reference=branch.get('baseline_branches');late_mask=phase_history[:count,0]>=warm+ramp+hold-.5
            if reference is None:
                branch['baseline_branches']=tuple(a.copy() for a in current);reference=branch['baseline_branches'];label='baseline'
            else:label='final' if record else None
            common=min(count,len(reference[0]));window=phase_history[:common,0]>=warm+ramp+hold-.5
            tests=dict(contact_topology=current[0][:common,15]!=reference[0][:common,15],wrapping_topology=np.any(current[2][:common]!=reference[2][:common],axis=1),muscle_branch=np.any(current[1][:common]!=reference[1][:common],axis=1))
            metrics['native_branch_comparison']=dict(first_late_difference_s={name:float(phase_history[np.flatnonzero(mask&window)[0],0]) if (mask&window).any() else None for name,mask in tests.items()},late_changed_frames={name:int(np.sum(mask&window)) for name,mask in tests.items()},maximum_late_toes_normal_difference_N=float(np.max(abs(current[0][:common,[8,10]][window]-reference[0][:common,[8,10]][window]))) if window.any() else None,maximum_late_MTP_velocity_difference_rad_s=float(np.max(abs(current[0][:common,3:5][window]-reference[0][:common,3:5][window]))) if window.any() else None,most_changed_muscles=[names[j] for j in np.argsort(np.sum(current[1][:common][window]!=reference[1][:common][window],axis=0))[-5:][::-1]] if window.any() else [])
            if label is not None:
                side='right' if branch is branches[0] else 'left';np.savez_compressed(out/(side+'_'+label+'_native_branches.npz'),phase=current[0],muscle_branch=current[1],wrap_counts=current[2],muscles=start['muscles'],phase_columns=np.array(['time','mtp_r','mtp_l','vmtp_r','vmtp_l','amtp_r','amtp_l','calcn_r_N','toes_r_N','calcn_l_N','toes_l_N','calcn_r_contacts','toes_r_contacts','calcn_l_contacts','toes_l_contacts','contact_crc32','solver_iter','nefc','wrap_crc32']))
        return residual,metrics,state,np.array(trace),point_residual
    def evaluate(z,force=False,record=False):
        nonlocal calls
        key=np.asarray(z).tobytes()
        if not force and not record and key in cache:return cache[key][0]
        if calls>=budget-1 and not record:raise RuntimeError('Declared paired rollout budget exhausted; final replay reserved')
        calls+=1;coefficients=np.array([modes@z[:nmode],modes@z[nmode:]]) if derivative_free_task else modes@z
        if full_family_audit or derivative_free_task:
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(max_workers=2) as executor:values=list(executor.map(lambda b:branch(b,coefficients,record),branches))
        else:values=[branch(b,coefficients,record) for b in branches]
        residual=np.r_[values[0][0],values[1][0],.01*z];cost=float(np.dot(residual,residual))
        if derivative_free_task:
            ratios=[]
            thresholds=dict(late_normal_error_BW=.002,final_dynamic_force_static_residual_max=.05,final_static_residual_max=.05,
                late_joint_speed_deg_s=.05,late_root_speed_mps=.0001,late_root_angular_speed_deg_s=.05,late_fiber_speed_l0_s=.001,
                late_force_variation_F0=1e-4,final_joint_error_deg=.25,final_root_error_m=.0005,final_root_rotation_error_deg=.25,
                terminal_hold_root_drift_m=.0005,terminal_hold_joint_drift_deg=.25,terminal_hold_root_rotation_drift_deg=.25,
                final_activation_excitation_error=1e-12,maximum_loaded_slip_mps=.005,maximum_root_translation_m=.005,
                maximum_root_rotation_deg=2.,source_acceleration_ratio=1.,source_velocity_ratio=1.,maximum_friction_cone_violation_N=1e-6)
            for value in values:
                metric=value[1];ratios.extend([metric[k]/v for k,v in thresholds.items()]);ratios.append(max(0.,-min(metric['minimum_mesh_gap_m']))/.002)
            ratios=np.nan_to_num(ratios,nan=1e6,posinf=1e6);terms=np.log1p(np.minimum(ratios,1e6)**2)
            cost=float(terms.max()+terms.mean()+sum(1e4*(2-v[1]['time_s']/(warm+ramp+hold)) for v in values if v[1]['guard']))
        history.append(dict(evaluation=calls,z=z.tolist(),shared_family_coefficients=coefficients.tolist(),cost=cost,branches=[v[1] for v in values]));cache[key]=(residual,values)
        if cost<best[0]:best[:]=[cost,z.copy(),values]
        (out/'SEARCH_HISTORY.json').write_text(json.dumps(history,indent=2))
        return residual
    candidate=json.loads(replay_candidate.read_text()) if replay_candidate is not None else None
    zero=np.array(candidate['z'],float) if candidate is not None else np.zeros(nmode*(2 if derivative_free_task else 1))
    if zero.shape!=(nmode*(2 if derivative_free_task else 1),) or not np.isfinite(zero).all() or max(abs(zero))>.1+1e-12:raise ValueError('Candidate outside declared shared-mode trust region')
    if candidate is not None:
        prior_basis=np.load(Path(candidate['source_audit'])/'shared_basis.npz')
        if np.max(abs(prior_basis['modes']-modes))>1e-12:raise ValueError('Candidate source basis changed')
    baseline=evaluate(zero).copy();first=cache[zero.tobytes()][1];repeat=evaluate(zero,force=True);parity=float(max(np.max(abs(a[2]-b[2])) for a,b in zip(first,cache[zero.tobytes()][1])))
    derivative=None;solver=None;reason=None
    if parity>1e-14:reason='Native augmented reset failed exact replay'
    elif derivative_free_task:
        def project(z):
            family=np.array([modes@z[:nmode],modes@z[nmode:]])
            return z/max(1.,float(max(abs(z)))/.1,float(max(np.sum(abs(family),axis=0)))/.1)
        proposal=np.zeros_like(zero)
        if proposal_audit is not None:
            prior=np.load(proposal_audit/'shooting_derivatives.npz');basis=np.load(proposal_audit/'shared_basis.npz')
            if not np.array_equal(basis['families'],families) or not np.array_equal(basis['phase_family_basis'],np.array([b['basis'] for b in branches])):raise ValueError('Proposal source-family basis changed')
            smooth=np.r_[0:23,46:138,584:607,630:722]
            smooth_J=prior['J'][smooth]@basis['modes'].T@modes
            proposal[:nmode]=np.linalg.lstsq(smooth_J,-baseline[smooth],rcond=1e-8)[0];proposal=project(proposal)
            manifest['proposal_scope']='Previously validated mean pose/fiber rows only; force/velocity gradients excluded. Every proposal scored by actual nonlinear dynamics.'
        original_cost=best[0]
        for scale in (.25,.5,1.):
            if np.linalg.norm(proposal)>1e-12:evaluate(project(scale*proposal))
        # Bounded task-coordinate pattern search; distal directions are mandatory.
        priorities=[3,0,4,6,5,1,2,7]
        for position,j in enumerate(priorities):
            anchor=best[1].copy()
            for sign in (1.,-1.):
                point=anchor.copy();point[j]+=sign*.05;evaluate(project(point))
            if position==3 and best[0]>=.995*original_cost:
                reason='No >=0.5% improvement after smooth proposals and hip/knee/ankle/MTP pattern checks';break
        if reason is None:
            direction=best[1][:nmode].copy()
            if max(abs(direction))<1e-12:direction=proposal[:nmode].copy()
            if max(abs(direction))>1e-12:
                direction*=.05/max(abs(direction));anchor=best[1].copy()
                for sign in (1.,-1.):
                    point=anchor.copy();point[nmode:]+=sign*direction;evaluate(project(point))
            reason='Completed bounded anatomical task-coordinate nonlinear pattern search; no physical-gradient fit'
        derivative=dict(method='Actual nonlinear derivative-free task search',smooth_proposal_only=proposal_audit is not None)
    elif candidate is not None:reason='Predeclared bounded nonlinear candidate check; no full-Jacobian fit'
    else:
        h=fd_step;J=np.column_stack([(evaluate(zero+np.eye(nmode)[j]*h)-baseline)/h for j in range(nmode)]);direction=np.arange(1.,nmode+1.);direction/=np.linalg.norm(direction)
        measured=evaluate(h*direction)-baseline;predicted=h*(J@direction);error=float(np.linalg.norm(measured-predicted)/max(1e-12,np.linalg.norm(measured)))
        derivative=dict(step=h,directional_relative_error=error,singular_values=np.linalg.svd(J,compute_uv=False).tolist())
        point=lambda z:np.r_[cache[np.asarray(z).tobytes()][1][0][4],cache[np.asarray(z).tobytes()][1][1][4],.01*z]
        point_base=point(zero);point_J=np.column_stack([(point(zero+np.eye(nmode)[j]*h)-point_base)/h for j in range(nmode)]);point_measured=point(h*direction)-point_base;point_predicted=h*(point_J@direction)
        np.savez_compressed(out/'shooting_derivatives.npz',J=J,baseline=baseline,directional_measured=measured,directional_predicted=predicted,point_J=point_J,point_measured=point_measured,point_predicted=point_predicted)
        derivative['point_directional_relative_error']=float(np.linalg.norm(point_measured-point_predicted)/max(1e-12,np.linalg.norm(point_measured)))
        if full_family_audit:
            half=evaluate(.5*h*direction)-baseline;double=evaluate(2*h*direction)-baseline;negative=evaluate(-h*direction)-baseline
            derivative['scale_checks']=dict(half_relative_error=float(np.linalg.norm(2*half-measured)/max(1e-12,np.linalg.norm(measured))),double_relative_error=float(np.linalg.norm(.5*double-measured)/max(1e-12,np.linalg.norm(measured))),opposite_relative_error=float(np.linalg.norm(-negative-measured)/max(1e-12,np.linalg.norm(measured))))
            np.savez_compressed(out/'shooting_scale_checks.npz',direction=direction,step=h,half=half,nominal=measured,double=double,negative=negative)
            reason='Full family derivative audit only; no fit authorized in this bounded job'
        elif error>.05:reason='Augmented shooting directional test failed; no derivative-based fit'
        else:
            try:
                def jacobian(z):
                    base=evaluate(z)
                    steps=np.where(z+h>.1,-h,h)
                    return np.column_stack([(evaluate(z+np.eye(nmode)[j]*steps[j])-base)/steps[j] for j in range(nmode)])
                solver=least_squares(evaluate,zero,jac=jacobian,bounds=(-.1,.1),max_nfev=6,ftol=1e-5,xtol=1e-5,gtol=1e-5)
            except RuntimeError as exception:reason=str(exception)
    winner=best[1];evaluate(winner,record=True);final=cache[winner.tobytes()][1]
    for side,value in zip(('right','left'),final):np.savez_compressed(out/(side+'_final_sparse.npz'),samples=value[3],nq=branches[0]['runner'].model.nq,nu=len(names),muscles=start['muscles'])
    final_coefficients=np.array([modes@winner[:nmode],modes@winner[nmode:]]) if derivative_free_task else modes@winner
    result=dict(accepted_load_bearing_transition=all(v[1]['accepted'] for v in final),accepted_dynamic_teacher=False,manifest=manifest,reset_parity_max=parity,derivative_audit=derivative,termination=reason,paired_rollout_evaluations=calls,wall_s=perf_counter()-started,total_native_steps=sum(b['runner'].physics_steps for b in branches),best_z=winner.tolist(),shared_family_coefficients=final_coefficients.tolist(),final=[v[1] for v in final])
    if candidate is not None:
        result['predeclared_candidate']=candidate
        result['actual_smooth_state_residual_norm']=float(np.linalg.norm(cache[winner.tobytes()][0][np.r_[0:23,46:138,584:607,630:722]]))
    (out/'MANIFOLD_SHOOTING.json').write_text(json.dumps(result,indent=2));np.savez_compressed(out/'shared_basis.npz',families=np.array(families),modes=modes,phase_family_basis=np.array([b['basis'] for b in branches]),source_torque_maps=np.array(torque_maps))
    print(json.dumps(result,indent=2))


def audit_augmented_inverse(model_folder, snapshot_file, reference_jacobian, spec_path, out, branch_derivatives=False):
    """Derivative validation,authority analysis,and one bounded command correction."""
    from time import perf_counter
    from scipy.optimize import lsq_linear
    from source_force_comparison import coordinate_stop_force
    if (out/'AUGMENTED_INVERSE_AUDIT.json').exists():
        raise ValueError('Preserve existing augmented inverse audit')
    with np.load(snapshot_file) as archive:
        z = {k: archive[k] for k in archive.files}
    with np.load(reference_jacobian) as archive:
        reference = {k: archive[k] for k in archive.files}
    m = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
    m.opt.timestep = float(z['timestep_s'])
    runner = AugmentedTransition(m, z['parameters'], float(z['filter_tau_s']))
    ar = CoupledArms(runner.model)
    selected = np.r_[np.arange(6), ar.rows]
    spec = json.loads(spec_path.read_text())
    addresses = np.array([v['qpos_address'] for v in spec['joints'].values()])
    low = np.array([v['native_bounds_rad'][0] for v in spec['joints'].values()])
    high = np.array([v['native_bounds_rad'][1] for v in spec['joints'].values()])
    runner.source_rom = (addresses, low, high)
    steps = int(z['window_steps']); results = []; arrays = []
    started = perf_counter()
    if not np.allclose(z['times'], reference['times'], atol=1e-12, rtol=0.):
        raise ValueError('Reference Jacobian times differ')
    for index, time in enumerate(z['times']):
        runner.load_snapshot(z, index)
        command = z['command'][index].copy()
        kp, kd = z['kp'][index], z['kd'][index]
        if not np.array_equal(kp[:7], kp[7:14]) or not np.array_equal(kd[:7], kd[7:14]):
            raise ValueError('Bilateral objective gains differ')
        qref, vref, aref = (z[key][index] for key in ('reference_q', 'reference_v', 'reference_a'))
        d = runner.data
        ar.update(d)
        E = np.zeros((m.nv, len(selected))); E[:6, :6] = np.eye(6); E[:, 6:] = ar.E
        mass = np.empty((m.nv, m.nv)); mujoco.mj_fullM(runner.model, mass, d.qM)
        reduced = E.T@mass@E
        root_weight = np.diag([abs(qref[1])]*3+[1.]*3)@reduced[:6, :6]

        def objective():
            model, data = runner.model, runner.data
            ar.update(data)
            constraint = data.qfrc_constraint
            if not bool(z['compensate_coordinate_stops']):
                constraint = constraint-coordinate_stop_force(model, data)
            target = ar.E.T@(data.qfrc_bias-data.qfrc_passive-constraint)+kp*(qref[ar.addresses]-data.qpos[ar.addresses])+kd*(vref[ar.rows]-data.qvel[ar.rows])
            actual = ar.arms@runner.backend.bank.force
            root_target = aref[:6]+25.*(qref[:6]-data.qpos[:6])+10.*(vref[:6]-data.qvel[:6])
            return np.r_[actual-target, root_weight@(data.qacc[:6]-root_target)]

        replay = runner.rollout(command, steps, z['window_excitation'][index])
        replay_error = float(np.max(abs(runner.data.qpos-z['window_qpos'][index, -1])))
        if replay_error > 1e-10:
            raise ValueError('Preallocated state restoration changed exact replay')
        baseline_state = runner.rollout(command, steps).copy()
        baseline = objective()
        baseline_activation_distance = runner.activation_distance[:steps].copy()
        baseline_fiber_velocity = runner.fiber_velocity[:steps].copy()
        state_J = np.zeros((len(baseline_state), m.nu)); acceleration_J = np.zeros((len(selected), m.nu))
        residual_J = np.zeros((len(baseline), m.nu)); fine_J = np.zeros_like(residual_J)
        for j in range(m.nu):
            for epsilon in (.005, .0025):
                lo = command.copy(); hi = command.copy()
                lo[j] = max(.01, command[j]-epsilon); hi[j] = min(.98, command[j]+epsilon)
                lower_state = runner.rollout(lo, steps).copy(); lower_a = runner.data.qacc[selected].copy(); lower_r = objective()
                upper_state = runner.rollout(hi, steps).copy(); upper_a = runner.data.qacc[selected].copy(); upper_r = objective()
                denominator = hi[j]-lo[j]
                derivative = (upper_r-lower_r)/denominator
                if epsilon == .005:
                    state_J[:, j] = (upper_state-lower_state)/denominator
                    acceleration_J[:, j] = (upper_a-lower_a)/denominator
                    residual_J[:, j] = derivative
                else:
                    fine_J[:, j] = derivative
        parity = float(np.max(abs(acceleration_J-reference['full_acceleration_J'][index])))
        half_error = float(np.linalg.norm(fine_J-residual_J)/max(1e-30, np.linalg.norm(fine_J)))
        directional = []
        rng = np.random.default_rng(731)
        for _ in range(2):
            direction = rng.uniform(-1., 1., m.nu)
            direction[command < .0102] = abs(direction[command < .0102])
            direction[command > .9798] = -abs(direction[command > .9798])
            h = 1e-4
            runner.rollout(command+h*direction, steps)
            observed = (objective()-baseline)/h
            expected = fine_J@direction
            directional.append(float(np.linalg.norm(observed-expected)/max(1e-30, np.linalg.norm(observed))))
        singular = np.linalg.svd(fine_J, compute_uv=False)
        rank = int(np.sum(singular > singular[0]*1e-6))
        dlow = np.maximum(.01-command, -.02); dhigh = np.minimum(.98-command, .02)
        reach_low = np.minimum(fine_J*dlow, fine_J*dhigh).sum(axis=1)
        reach_high = np.maximum(fine_J*dlow, fine_J*dhigh).sum(axis=1)
        unreachable = np.flatnonzero((-baseline < reach_low) | (-baseline > reach_high))
        derivative_valid = parity <= 1e-6*(1.+np.max(abs(reference['full_acceleration_J'][index]))) and half_error <= .05 and max(directional) <= .05
        branch_report = None
        if branch_derivatives:
            # At zero smoothing the activation slope changes when excitation
            # crosses activation. Test directional derivatives, not a fictitious
            # smooth central Jacobian. This mode never solves an inverse.
            sides = np.zeros((2, len(baseline), m.nu))
            half_sides = np.zeros_like(sides)
            available = np.zeros((2, m.nu), dtype=bool)
            activation_switches = np.zeros(2, dtype=int)
            fiber_switches = np.zeros(2, dtype=int)
            for side, sign in enumerate((-1., 1.)):
                for j in range(m.nu):
                    room = command[j]-.01 if sign < 0 else .98-command[j]
                    if room < 1e-4:
                        continue
                    available[side, j] = True
                    for epsilon, destination in ((1e-4, sides), (5e-5, half_sides)):
                        perturbation = command.copy()
                        perturbation[j] += sign*epsilon
                        runner.rollout(perturbation, steps)
                        destination[side, :, j] = (objective()-baseline)/(sign*epsilon)
                        if epsilon == 1e-4:
                            activation_switches[side] += int(np.sum((runner.activation_distance[:steps]>0.) != (baseline_activation_distance>0.)))
                            fiber_switches[side] += int(np.sum((runner.fiber_velocity[:steps]>0.) != (baseline_fiber_velocity>0.)))
            branch_directional = []
            central_small_directional = []
            rng = np.random.default_rng(731)
            central_small = .5*(half_sides[0]+half_sides[1])
            central_small[:, ~available[0]] = half_sides[1, :, ~available[0]].T
            central_small[:, ~available[1]] = half_sides[0, :, ~available[1]].T
            for _ in range(2):
                direction = rng.uniform(-1., 1., m.nu)
                direction[~available[0]] = abs(direction[~available[0]])
                direction[~available[1]] = -abs(direction[~available[1]])
                h = 1e-5
                runner.rollout(command+h*direction, steps)
                observed = (objective()-baseline)/h
                selected_columns = np.where(direction[None, :]>=0., half_sides[1], half_sides[0])
                expected = selected_columns@direction
                denominator = max(1e-30, np.linalg.norm(observed))
                branch_directional.append(float(np.linalg.norm(observed-expected)/denominator))
                central_small_directional.append(float(np.linalg.norm(observed-central_small@direction)/denominator))
            branch_report = dict(one_sided_step=1e-4, half_step=5e-5, directional_step=1e-5,
                                 half_step_relative_error=float(np.linalg.norm(sides-half_sides)/max(1e-30, np.linalg.norm(half_sides))),
                                 signed_column_directional_errors=branch_directional,
                                 small_central_directional_errors=central_small_directional,
                                 activation_switch_counts=activation_switches.tolist(), fiber_velocity_sign_switch_counts=fiber_switches.tolist(),
                                 nominal_activation_distance_below_1e8=int(np.sum(abs(baseline_activation_distance)<1e-8)),
                                 nominal_fiber_velocity_below_1e8=int(np.sum(abs(baseline_fiber_velocity)<1e-8)),
                                 scope='Derivative diagnosis only. No source-law smoothing, inverse, or physical-controller trial.')
        inverse = dict(attempted=False)
        if derivative_valid and not branch_derivatives:
            regularization = .05*singular[0]
            fit = lsq_linear(np.vstack([fine_J, regularization*np.eye(m.nu)]), np.r_[-baseline, np.zeros(m.nu)],
                             bounds=(dlow, dhigh), method='bvls', tol=1e-10, max_iter=200)
            delta = fit.x
            runner.rollout(command+delta, steps)
            actual = objective()
            excess = max(float(np.max(low-runner.data.qpos[addresses])), float(np.max(runner.data.qpos[addresses]-high)), 0.)
            inverse = dict(attempted=True, regularization=regularization, command_correction_limit=.02,
                           maximum_command_change=float(np.max(abs(delta))), original_objective_norm=float(np.linalg.norm(baseline)),
                           predicted_objective_norm=float(np.linalg.norm(baseline+fine_J@delta)), actual_objective_norm=float(np.linalg.norm(actual)),
                           nonlinear_objective_reduced=bool(np.linalg.norm(actual)<np.linalg.norm(baseline)), source_ROM_excess_rad=excess,
                           scope='One bounded regularized correction,not parameter optimization or a physical teacher acceptance.')
        results.append(dict(time_s=float(time), exact_preallocated_replay_qpos_error=replay_error,
                            column_parity_max_error=parity, half_step_J_relative_error=half_error,
                            independent_directional_relative_errors=directional, derivative_gates_pass=bool(derivative_valid),
                            singular_values=singular.tolist(), effective_rank_1ppm=rank,
                            condition_number=float(singular[0]/max(1e-30, singular[-1])),
                            unreachable_rows_with_small_box=unreachable.tolist(), inverse=inverse,
                            one_sided_branch_diagnostics=branch_report,
                            residual_rows=list(JOINTS)+['root_'+str(i) for i in range(6)],
                            baseline_residual=baseline.tolist(), reachable_correction_low=reach_low.tolist(), reachable_correction_high=reach_high.tolist()))
        arrays.append((state_J, residual_J, fine_J, baseline))
    result = dict(source_snapshots=str(snapshot_file), reference_jacobian=str(reference_jacobian), horizon_s=steps*float(z['timestep_s']),
                  rows=results, wall_s=perf_counter()-started, native_physics_steps=runner.physics_steps,
                  objective='Existing shared bilateral joint PD torque and root5rad/s feedback,evaluated on complete future native state.',
                  scope='No physical model,reference endpoint,gain,or force floor changes. Derivatives validated before a small regularized source-command correction; no live teacher or Hunt C fitting.')
    out.mkdir(parents=True, exist_ok=True)
    (out/'AUGMENTED_INVERSE_AUDIT.json').write_text(json.dumps(result, indent=2))
    np.savez_compressed(out/'augmented_derivatives.npz', times=z['times'], augmented_state_J=np.array([r[0] for r in arrays]),
                        objective_J=np.array([r[1] for r in arrays]), fine_objective_J=np.array([r[2] for r in arrays]), baseline_residual=np.array([r[3] for r in arrays]))
    print(json.dumps(result, indent=2))


def audit_full_state_jacobian(model_folder, snapshot_file, out):
    """Compare conditional force inversion with an independent full-body rollout."""
    import copy
    from thelen_dynamic import ConstitutiveBackend
    from source_force_comparison import NativeAccelerationTeacher
    if (out/'FULL_ACTUATION_JACOBIAN.json').exists():
        raise ValueError('Preserve existing full-state audit')
    with np.load(snapshot_file) as archive:
        z = {key: archive[key] for key in archive.files}
    model = mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
    model.opt.timestep = float(z['timestep_s'])
    mask = mujoco.mjtState(int(z['state_spec']))
    dt = float(z['timestep_s']); steps = int(z['window_steps'])
    parameters = z['parameters']; beta = 1.-np.exp(-dt/float(z['filter_tau_s']))
    rows = []; matrices = []

    def restore(index):
        m = copy.copy(model); d = mujoco.MjData(m)
        mujoco.mj_setState(m, d, z['integration_state'][index], mask)
        mujoco.mj_forward(m, d)
        backend = ConstitutiveBackend(m, d, parameters)
        backend.bank.state[:] = z['fiber_state'][index]
        backend.bank.force[:] = z['source_force'][index]
        m.actuator_biasprm[:] = z['actuator_biasprm'][index]
        mujoco.mj_setState(m, d, z['integration_state'][index], mask)
        mujoco.mj_forward(m, d)
        return m, d, backend

    def rollout(index, command, recorded=False):
        m, d, backend = restore(index)
        excitation = z['excitation'][index].copy()
        for k in range(steps):
            if recorded:
                d.ctrl[:] = z['window_excitation'][index, k]
            else:
                excitation += beta*(command-excitation)
                d.ctrl[:] = excitation
            backend.step()
        return d.qpos.copy(), d.qvel.copy(), d.act.copy(), backend.bank.force.copy(), d.qacc.copy()

    for index, time in enumerate(z['times']):
        m, d, backend = restore(index)
        ar = CoupledArms(m); ar.update(d)
        selected = np.r_[np.arange(6), ar.rows]
        clone_error = float(np.max(abs(d.qacc-z['initial_qacc'][index])))
        replay = rollout(index, z['command'][index], recorded=True)
        keys = ('qpos', 'qvel', 'activation', 'force', 'qacc')
        replay_error = {key: float(np.max(abs(replay[k]-z['window_'+key][index, -1]))) for k, key in enumerate(keys)}
        if clone_error > 1e-7 or replay_error['qpos'] > 1e-10 or replay_error['qvel'] > 1e-9 or replay_error['force'] > 1e-6:
            raise ValueError('Full augmented-state restore/replay failed: '+str(replay_error))
        current = d.ctrl.copy(); current_acceleration = d.qacc.copy()
        d.ctrl[:] = np.minimum(.98, current+.005); mujoco.mj_forward(m, d)
        direct_feedthrough = float(np.max(abs(d.qacc-current_acceleration)))
        d.ctrl[:] = current; mujoco.mj_forward(m, d)
        predictor = SourcePredictor(parameters, m.actuator_dynprm, dt, float(z['filter_tau_s']))
        forecast = predictor.linearize_command(backend.bank.state, d.act, z['excitation'][index], z['command'][index], d.actuator_length, d.actuator_velocity, steps)
        if forecast is None:
            raise ValueError('Conditional forecast failed at audit state')
        predicted_force, gain, command = forecast
        mapper = NativeAccelerationTeacher(m, ar, np.ones(m.nu, bool))
        mapper.update(d, gain, d.act, d.qpos.copy(), d.qvel.copy(), np.zeros(m.nv), command,
                      source_prediction=(predicted_force, command, backend.bank.force.copy()))
        approximate = mapper.G*gain[None, :]
        full = np.zeros_like(approximate); force_jacobian = np.zeros((m.nu, m.nu))
        for j in range(m.nu):
            lo = command.copy(); hi = command.copy()
            lo[j] = max(.01, command[j]-.005); hi[j] = min(.98, command[j]+.005)
            lower = rollout(index, lo); upper = rollout(index, hi)
            denominator = hi[j]-lo[j]
            full[:, j] = (upper[4][selected]-lower[4][selected])/denominator
            force_jacobian[:, j] = (upper[3]-lower[3])/denominator
        mass = np.empty((m.nv, m.nv)); mujoco.mj_fullM(m, mass, d.qM)
        E = np.zeros((m.nv, len(selected))); E[:6, :6] = np.eye(6); E[:, 6:] = ar.E
        metric = np.linalg.cholesky(E.T@mass@E).T
        weighted_full = metric@full; weighted_approximate = metric@approximate
        column_norm = np.linalg.norm(weighted_full, axis=0)
        significant = column_norm > 1e-6*column_norm.max()
        cosine = np.sum(weighted_full*weighted_approximate, axis=0)/np.maximum(1e-30, column_norm*np.linalg.norm(weighted_approximate, axis=0))
        off_diagonal = force_jacobian-np.diag(np.diag(force_jacobian))
        rows.append(dict(time_s=float(time), clone_acceleration_max_error=clone_error, exact_window_replay_errors=replay_error,
                         current_acceleration_command_feedthrough_max=direct_feedthrough,
                         weighted_acceleration_J_relative_error=float(np.linalg.norm(weighted_approximate-weighted_full)/np.linalg.norm(weighted_full)),
                         significant_columns=int(significant.sum()), opposite_direction_columns=int(np.sum(cosine[significant]<0.)),
                         force_J_offdiagonal_relative_norm=float(np.linalg.norm(off_diagonal)/np.linalg.norm(force_jacobian)),
                         force_J_diagonal_relative_error=float(np.linalg.norm(np.diag(force_jacobian)-gain)/np.linalg.norm(force_jacobian)),
                         native_map_validation=mapper.diagnostics()))
        matrices.append((full, approximate, force_jacobian, gain))
    result = dict(source_snapshots=str(snapshot_file), timestep_s=dt, horizon_s=steps*dt, rows=rows,
                  scope='Independent full-native finite differences propagate q,v,activation,fibers,and excitation filtering. Conditional forecast holds only a linear length prediction. No controller/physics parameters are fitted or changed.')
    out.mkdir(parents=True, exist_ok=True)
    (out/'FULL_ACTUATION_JACOBIAN.json').write_text(json.dumps(result, indent=2))
    np.savez_compressed(out/'jacobian_windows.npz', times=z['times'], full_acceleration_J=np.array([r[0] for r in matrices]),
                        approximate_acceleration_J=np.array([r[1] for r in matrices]), full_force_J=np.array([r[2] for r in matrices]), conditional_force_gain=np.array([r[3] for r in matrices]))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('model', 'out'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--trace', type=Path)
    p.add_argument('--snapshots', type=Path)
    p.add_argument('--inverse-reference-jacobian', type=Path)
    p.add_argument('--spec', type=Path)
    p.add_argument('--branch-derivatives', action='store_true')
    p.add_argument('--stationary-start', type=Path)
    p.add_argument('--stationary-target', type=Path)
    p.add_argument('--verify-stationary-inverse', type=Path)
    p.add_argument('--physics-dt', type=float)
    p.add_argument('--continuation-seed', type=Path)
    p.add_argument('--verify-continuation', type=Path)
    p.add_argument('--audit-transfer-geometry', type=Path)
    p.add_argument('--rolling-stance', action='store_true')
    p.add_argument('--validate-contact-kinematics', type=Path, nargs='+')
    p.add_argument('--extend-transfer-from', type=Path)
    p.add_argument('--audit-transfer-authority', type=Path)
    p.add_argument('--prepare-stationary-transfer', type=Path)
    p.add_argument('--stationary-source-capacity', type=Path)
    p.add_argument('--stationary-manifold-seed',type=Path)
    p.add_argument('--load-fraction',type=float,default=.01)
    p.add_argument('--shoot-stationary-manifold',type=Path,nargs=2)
    p.add_argument('--paired-rollout-budget',type=int,default=24)
    p.add_argument('--terminal-window',action='store_true')
    p.add_argument('--shooting-candidate',type=Path)
    p.add_argument('--full-family-audit',action='store_true')
    p.add_argument('--derivative-free-task',action='store_true')
    p.add_argument('--proposal-audit',type=Path)
    p.add_argument('--finite-difference-step',type=float,default=.0001)
    p.add_argument('--samples', type=int, default=101)
    a = p.parse_args()
    if a.rolling_stance and a.continuation_seed is None:
        p.error('The opt-in rolling stance task requires --continuation-seed')
    if a.extend_transfer_from is not None and not a.rolling_stance:
        p.error('Extended transfer requires the explicitly selected rolling-stance task')
    if a.shoot_stationary_manifold:
        if a.stationary_start is None or a.spec is None:p.error('Shooting requires stationary start and frozen spec')
        shoot_stationary_manifold(a.model,a.stationary_start,a.shoot_stationary_manifold,a.spec,a.out,a.paired_rollout_budget,a.physics_dt,a.terminal_window,a.shooting_candidate,a.full_family_audit,a.finite_difference_step,a.derivative_free_task,a.proposal_audit)
    elif a.stationary_manifold_seed:
        if a.stationary_start is None:p.error('Manifold seed requires --stationary-start')
        prepare_manifold_seed(a.model,a.stationary_start,a.stationary_manifold_seed,a.out,a.load_fraction)
    elif a.stationary_source_capacity:
        stationary_source_capacity(a.model,a.stationary_source_capacity,a.out)
    elif a.prepare_stationary_transfer:
        if a.stationary_start is None: p.error('Stationary preparation requires --stationary-start')
        prepare_stationary_transfer(a.model,a.stationary_start,a.prepare_stationary_transfer,a.out)
    elif a.audit_transfer_authority:
        if a.stationary_start is None: p.error('Authority audit requires --stationary-start')
        audit_transfer_authority(a.model,a.stationary_start,a.audit_transfer_authority,a.out)
    elif a.validate_contact_kinematics:
        validate_contact_kinematics(a.model,a.validate_contact_kinematics,a.out)
    elif a.audit_transfer_geometry:
        if a.stationary_start is None:
            p.error('Transfer geometry audit requires --stationary-start')
        audit_transfer_geometry(a.model,a.stationary_start,a.audit_transfer_geometry,a.out)
    elif a.stationary_start:
        if a.stationary_target is None or a.spec is None:
            p.error('Stationary group inverse requires --stationary-target and --spec')
        stationary_group_inverse(a.model, a.stationary_start, a.stationary_target, a.spec, a.out,a.verify_stationary_inverse,a.physics_dt,a.continuation_seed,a.verify_continuation,a.rolling_stance,a.extend_transfer_from)
    elif a.snapshots and a.inverse_reference_jacobian:
        if a.spec is None:
            p.error('Augmented inverse audit requires --spec')
        audit_augmented_inverse(a.model, a.snapshots, a.inverse_reference_jacobian, a.spec, a.out, a.branch_derivatives)
    elif a.snapshots:
        audit_full_state_jacobian(a.model, a.snapshots, a.out)
    elif a.trace:
        audit(a.model, a.trace, a.out, a.samples)
    else:
        p.error('Provide --trace or --snapshots')
