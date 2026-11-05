"""Persistent Rustuna worker for ordered *offline neural* Hunt A-D fits.

Reference joint/muscle samples are replayed, not mechanically integrated. No
joint tracking, foot clearance, physical stage pass, or walking is inferred.
Import is inert. First real preflight is deliberately deferred until Ben frees
the machine. All optimized parameters are shared across left/right sides.
"""
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.signal import find_peaks

S = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(S))
from controller_profiles_20261104 import profile as P
from controller_profiles_20261104.prepare import audit
from muscle_map import classify

INSTANCE = None


class Replay:
    def __init__(self, manifest):
        self.manifest = manifest
        self.dt = .002
        path = Path(manifest['targets'])
        if hashlib.sha256(path.read_bytes()).hexdigest() != manifest['targets_sha256']:
            raise ValueError('Reference targets changed')
        z = np.load(path, allow_pickle=False)
        self.acts = z['act_names'].tolist()
        self.net = P.build(manifest['family'], self.acts, options=manifest.get('hunt_options'), config=manifest['controller_profile'])
        self.static = audit(self.net)
        if not self.static['static_pass']:
            raise ValueError('Static wiring gate failed')
        self.phase = manifest['phase']
        self.T = float(manifest['target_period_s'])
        duration = (8. * self.T if self.phase == 'A' else 12. * self.T)
        self.t = np.arange(round(duration / self.dt)) * self.dt
        query = manifest['target_start_s'] + self.t % self.T
        if query.min() < z['t'].min() or query.max() > z['t'].max():
            raise ValueError('Reference interval outside target data')
        def interp(key):
            a = z[key]
            return np.stack([np.interp(query, z['t'], a[:, j]) for j in range(a.shape[1])], axis=1)
        self.desired = interp('activation')
        length, velocity, force = interp('length'), interp('velocity'), interp('force')
        mid = np.median(z['length'], axis=0)
        half = np.maximum(np.ptp(z['length'], axis=0) / 2., .02)
        import params
        self.u = np.zeros((len(self.t), len(self.net.inputs)))
        ports = self.net.command_runtime.ports
        # Offline expected-signal adapter: declared physical scales, not an
        # assertion that these reference forces are measured animal afferents.
        encoded = dict(Ia=params.AFF['ia_gain'] * np.clip(velocity / .6, -1., 1.),
                       II=params.AFF['ii_gain'] * np.clip((length - mid) / half, -1., 1.),
                       Ib=params.AFF['ib_gain'] * np.maximum(force, 0.) / np.maximum(z['Fmax'], 1.))
        for i, act in enumerate(self.acts):
            for aff in ('Ia', 'II', 'Ib'):
                if aff + '_' + act in ports:
                    self.u[:, ports[aff + '_' + act]] = encoded[aff][:, i]
        for side in ('l', 'r'):
            selected = [i for i, act in enumerate(self.acts) if act.endswith('_' + side) and classify(act).groups[0] in ('hip_ext', 'knee_ext', 'ankle_pf')]
            self.u[:, ports['LOAD_c_' + side]] = np.clip(np.sum(np.maximum(force[:, selected], 0.), axis=1) / max(float(np.percentile(np.sum(np.maximum(force[:, selected], 0.), axis=1), 95)), 1.), 0., 1.)
        self.c = self.net.compiled
        self.core = self.c.c
        if hasattr(self.net, 'apply_weights'):
            self.net.apply_weights(1.)
        self.initial = {k: getattr(self.core, k).copy() for k in ('V', 'b_gate', 'c_gate', 'tau_max_b', 'g_ion')}
        self.base_g = self.c.gmax.copy()
        self.names = P.nodes(self.net)
        self.rg = np.array([self.net.idx['RG_' + phase + '_' + side] for side in ('l', 'r') for phase in ('E', 'F')])
        self.mn = np.array([self.net.idx['MN_' + act] for act in self.acts])
        self.pf = np.array([i for i, n in enumerate(self.names) if n.startswith('PF_') and not n.startswith(('PF_IN', 'PF_HUNT_TONIC'))])
        self.selected = np.array([classify(act).groups[0] in ('hip_ext', 'hip_flex', 'knee_ext', 'knee_flex', 'ankle_pf', 'ankle_df') for act in self.acts])
        self.output = np.empty((len(self.t), len(self.acts)))
        self.rg_trace = np.empty((len(self.t), 4))
        self.pf_range_min = np.empty(len(self.pf)); self.pf_range_max = np.empty(len(self.pf))
        self.pf_trace = np.empty((len(self.t), len(self.pf)))
        self.work = np.zeros(len(self.net.inputs))

    def evaluate(self, parameters, deadline=None):
        tick = time.perf_counter()
        self.c.gmax[:] = self.base_g
        for key, value in self.initial.items():
            getattr(self.core, key)[:] = value
        self.core.tau_max_b[:, self.rg] *= parameters.get('tau_scale', 1.)
        self.core.g_ion[:, self.rg] *= parameters.get('g_ion_scale', 1.)
        self.c.refresh_parameters()
        runtime = self.net.command_runtime
        runtime.unloaded.fill(0.); runtime.armed.fill(False)
        runtime.loaded.fill(False); runtime.heel_trace.fill(0.)
        runtime.sensory_enabled = self.phase != 'A'
        runtime.contact_enabled = False
        for key in runtime.shared_gains:
            runtime.shared_gains[key] = parameters.get(key, 1.)
        for i, (pre, post) in enumerate(zip(self.c.pre, self.c.post)):
            source, target = self.names[pre], self.names[post]
            if source.startswith('V3E_') or (source.startswith('RG_E_') and target.startswith('V3E_')):
                self.c.gmax[i] *= parameters.get('commissural_E',1.)
            elif source.startswith(('CIN_F_','V2a_','V0V_','V0D_','INI_')) or (source.startswith('RG_F_') and target.startswith(('CIN_F_','V2a_','V0D_'))):
                self.c.gmax[i] *= parameters.get('commissural_F',1.)
            if target.startswith('MN_') and source.startswith('PF_') and not source.startswith('PF_IN'):
                group = classify(target[3:]).groups[0]
                self.c.gmax[i] *= parameters.get('ff_' + group, 1.)
            if not P.central(target):
                for aff in ('Ia', 'II', 'Ib'):
                    if source.startswith(aff + '_'):
                        self.c.gmax[i] *= parameters.get('local_' + aff, 1.)
            if target.startswith('MN_'):
                if source.startswith(('IBIN_', 'IBEXC_')):
                    self.c.gmax[i] *= parameters.get('relay_Ib_output', 1.)
                elif source.startswith(('IIX_', 'IIIN_')):
                    self.c.gmax[i] *= parameters.get('relay_II_output', 1.)
        self.pf_range_min.fill(np.inf); self.pf_range_max.fill(-np.inf)
        for i, t in enumerate(self.t):
            if i % 128 == 0 and deadline is not None and time.time() >= deadline:
                raise TimeoutError('Neural replay deadline')
            np.copyto(self.work, self.u[i])
            # One initial command pulse, then genuinely undriven A. Later
            # blocks retain a declared symmetric physiological command level.
            for side in ('l', 'r'):
                runtime.commands['F_' + side] = 0. if self.phase == 'A' else self.manifest['command_F']
                runtime.commands['E_' + side] = 0. if self.phase == 'A' else self.manifest['command_E']
            if t < .1:
                runtime.commands['E_r'] = min(1.,runtime.commands['E_r']+.4)
                runtime.commands['F_l'] = min(1.,runtime.commands['F_l']+.4)
            runtime.prepare_inputs(self.work)
            try:
                self.c.forward(self.work)
            except (FloatingPointError, OverflowError):
                return dict(ok=True, value=1e6, diagnostics=dict(finite=False, accepted=False,
                    failure='nonfinite_neural_trial', stopped_step=i,
                    scope='Offline neural replay; not a body result'))
            self.rg_trace[i] = self.core.V[self.rg]
            self.pf_trace[i] = self.core.V[self.pf]
            self.output[i] = np.clip(self.core.V[self.mn] / 5., 0., 1.)
            if t >= 2. * self.T:
                np.minimum(self.pf_range_min, self.core.V[self.pf], out=self.pf_range_min)
                np.maximum(self.pf_range_max, self.core.V[self.pf], out=self.pf_range_max)
        keep = self.t >= 2. * self.T
        periods, cycles, peak_times = [], [], []
        for column in (1, 3):
            peaks, _ = find_peaks(self.rg_trace[keep, column], prominence=.15, distance=max(1, int(.3 * self.T / self.dt)))
            cycles.append(len(peaks))
            peak_times.append(peaks * self.dt)
            periods.append(float(np.median(np.diff(peaks)) * self.dt) if len(peaks) >= 3 else None)
        target = self.T * (2. if self.phase == 'A' else 1.)
        finite = bool(np.isfinite(self.output).all() and np.isfinite(self.rg_trace).all())
        period_error = float(np.mean([abs(p / target - 1.) if p else 10. for p in periods]))
        phase_error = 1.
        if all(p is not None for p in periods) and all(len(p) >= 3 for p in peak_times):
            separation = np.median(np.min(np.abs(peak_times[0][:, None] - peak_times[1][None, :]), axis=1)) / np.mean(periods)
            phase_error = float(abs(min(separation, .5) - .5))
        rmse = float(np.sqrt(np.mean((self.output[keep][:, self.selected] - self.desired[keep][:, self.selected]) ** 2)))
        peak_error = float(np.sqrt(np.mean((self.output[keep][:, self.selected].max(axis=0) - self.desired[keep][:, self.selected].max(axis=0)) ** 2)))
        pf_dynamic = int(np.sum(self.pf_range_max - self.pf_range_min > .05))
        # Zero/frozen motion is explicitly a failure, regardless of MN RMSE.
        rhythmic = finite and min(cycles) >= 3 and period_error <= .10 and phase_error <= .15
        value = (1e6 if not finite else 1e4 if min(cycles) < 3 else 0.) + 20. * period_error + 20. * phase_error
        if self.phase == 'C':
            value += 10. * peak_error
        elif self.phase == 'D':
            value += 10. * rmse
        diagnostics = dict(finite=finite, rg_cycles=cycles, rg_period_s=periods,
                           target_period_s=target, relative_period_error=period_error,
                           bilateral_phase_error=phase_error,
                           mn_command_rmse=rmse, dynamic_pf_cells=pf_dynamic,
                           mn_peak_rmse=peak_error,
                           accepted=bool(rhythmic and (self.phase not in ('C', 'D') or (peak_error if self.phase == 'C' else rmse) <= .25 and pf_dynamic > 0)),
                           wall_s=time.perf_counter() - tick,
                           scope='Offline neural replay on expected afferents; no physics integration, joint tracking or curriculum acceptance')
        return dict(ok=True, value=float(value), diagnostics=diagnostics)


def handle_request(request):
    global INSTANCE
    if request['op'] == 'preflight':
        INSTANCE = Replay(request['manifest'])
        baseline = INSTANCE.evaluate(request['manifest'].get('fixed_parameters', {}), request.get('deadline_epoch_s'))
        return dict(ok=True, ready_for_optimization=baseline['diagnostics']['finite'],
                    baseline=baseline, static=INSTANCE.static,
                    runtime_verification='First neural replay occurs here, not during preparation')
    if INSTANCE is None or request['op'] not in ('evaluate','capture'):
        raise ValueError('Preflight required')
    result=INSTANCE.evaluate(request['params'], request.get('deadline_epoch_s'))
    if request['op']=='capture':
        path=Path(request['path']).resolve()
        root=Path(r'D:\Users\Ben Bolen\ControllerTraining_20261104').resolve()
        if not path.is_relative_to(root):raise ValueError('Capture outside owned training root')
        np.savez_compressed(path,t=INSTANCE.t,rg=INSTANCE.rg_trace,pf=INSTANCE.pf_trace,
            mn=INSTANCE.output,desired=INSTANCE.desired,
            rg_names=np.array([INSTANCE.names[i]for i in INSTANCE.rg]),
            pf_names=np.array([INSTANCE.names[i]for i in INSTANCE.pf]),acts=np.array(INSTANCE.acts))
    return result
