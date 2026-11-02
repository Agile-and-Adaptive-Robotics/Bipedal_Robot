"""One distinct ideal-force native static capacity case; no coactivation claim."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear, minimize
import cad_activation_semantics as A
import cad_rigid_stance as C

OUT = A.OUT


class ZeroFloor(C.Native):
    def inverse_fit(self, x):
        heights = self.pose(x)
        model, data = self.m, self.d
        data.qacc[:] = 0
        C.mujoco.mj_inverse(model, data)
        rig = np.zeros(model.nv)
        C.mujoco.mj_applyFT(model, data, [0, .85*self.weight, 0], [0, 0, 0],
                          data.xipos[self.pelvis], self.pelvis, rig)
        matrix = -self.s.dynamic.length_jacobian[self.s.tids].T
        target = data.qfrc_inverse - rig
        fit = lsq_linear(self.scale[:, None] * matrix * self.s.caps,
                         self.scale * target, bounds=(0., 1.),
                         method='bvls', tol=1e-12, max_iter=200)
        self.force[:] = fit.x * self.s.caps
        error = self.scale * (matrix @ self.force + rig - data.qfrc_inverse)
        gaps = np.asarray(self.s.constraints())
        return error, np.maximum(gaps-.0005, 0)/.00005, np.maximum(heights+.00001, 0)/.0001, np.maximum(-.004-heights, 0)/.0001, heights, fit


def main():
    audit = json.loads((OUT/'ACTIVATION_SEMANTICS_AUDIT.json').read_text())
    assert not audit['physiological_5percent_activation_map_available']
    initial = C.OUT/'REFINED_PRELOAD_STATE.npz'
    protocol = dict(machine='easteregg2', scope='NEW ideal54 zero-floor mechanical static capacity case, not SNS or BPA or 5percent coactivation',
                    model_definition=audit['zero_floor_fixture_definition'],
                    previous_floor_negatives_preserved=True,
                    initial_state=str(initial), initial_sha256=hashlib.sha256(initial.read_bytes()).hexdigest(),
                    numerics=dict(solref=[.004, 1], impedance=.9984069577958117, timestep=.002, iterations=500),
                    geometry_mass_damping_rig_capacity_unchanged=True,
                    support=.85, source_floor_change='Explicit new mechanical bounds: all54 T in [0,F0], no passive force model.',
                    maxfev=1500, method='Powell, one bounded case with existing pose-family bounds',
                    gates=dict(qacc=.05, COMacc_norm=.01, ground_error_N=.1, each_foot_N=1., each_heel_toe_N=.5,
                               loop_gap_m=.0005, maximum_knee_rad=-.02, maximum_heading_deg=10., mesh_penetration_m=.004, EOM=1e-6),
                    no_physiological_or_neural_success_claim=True,
                    source_sha256={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in (Path(__file__),Path(A.__file__),Path(C.__file__))})
    path = OUT/'PREREGISTRATION.json'
    if path.exists():
        assert json.loads(path.read_text()) == protocol
    else:
        path.write_text(json.dumps(protocol, indent=2))
    b = ZeroFloor()
    calls = [0]
    def objective(x):
        error, gap, air, depth, _, _ = b.inverse_fit(x)
        calls[0] += 1
        residual = np.r_[error, gap, air, depth, 1e-8*x]
        return float(residual@residual)
    start = C.P.load(initial)['variables']
    lower = np.r_[np.full(6, -10.), .1, np.full(4, -4.)]
    upper = np.r_[np.full(6, 10.), 4., np.full(4, 4.)]
    clock = time.perf_counter()
    fit = minimize(objective, start, method='Powell', bounds=list(zip(lower, upper)),
                   options=dict(maxfev=1500, xtol=1e-9, ftol=1e-12))
    error, _, _, _, heights, _ = b.inverse_fit(fit.x)
    com_acc, loads = b.actual()
    qacc = b.d.qacc.copy()
    gaps = b.s.constraints()
    knees = list(b.s.dynamic.coordinates(b.d).values())
    headings = []
    for side in ('r', 'l'):
        axis = b.d.xmat[b.m.body('foot_'+side).id].reshape(3, 3)[:, 0]
        headings.append(float(np.rad2deg(np.arctan2(axis[2], axis[0]))))
    ledger = C.accounting(b)
    gates = dict(qacc=np.max(np.abs(qacc)) <= .05, COMacc=np.linalg.norm(com_acc) <= .01,
                 ground=abs(loads.sum()-.15*b.weight) <= .1,
                 both_feet=np.all(loads.sum(axis=1) > 1), heel_toe=np.all(loads > .5),
                 loop=max(gaps) <= .0005, knees=max(knees) <= -.02,
                 heading=max(np.abs(headings)) <= 10,
                 mesh=np.all(heights<=0) and np.all(heights>=-.004),
                 EOM=ledger['native_EOM_residual_max'] < 1e-6,
                 physical_rank4=ledger['physical_rank'] == 4)
    gates = {key:bool(value) for key,value in gates.items()}
    result = dict(accepted=all(gates.values()), gates=gates, qacc_max=float(np.max(np.abs(qacc))),
                  COM_acceleration=com_acc.tolist(), ground_N=float(loads.sum()), target_ground_N=.15*b.weight,
                  contact_loads_N=loads.tolist(), loop_gaps_m=gaps, knees_rad=knees, heading_deg=headings,
                  foot_min_y_m=heights.tolist(), force_N=b.force.tolist(), native_accounting=ledger,
                  scaled_inverse_residual_max=float(np.max(np.abs(error))), solver_success=bool(fit.success),
                  nfev=int(fit.nfev), inverse_calls=calls[0]+1, wall_s=time.perf_counter()-clock,
                  neural_or_pressure_or_physiological_coactivation=False,
                  case='Explicit ideal-force zero-floor capacity fixture', numerical_parameters_frozen=True)
    np.savez_compressed(OUT/'PRELOAD_STATE.npz', qpos=b.d.qpos.copy(), qvel=b.d.qvel.copy(),
                        qacc=qacc, force_N=b.force.copy(), variables=fit.x, root_anchor=b.d.qpos[:3].copy())
    (OUT/'STATIC_RESULT.json').write_text(json.dumps(result, indent=2))
    print(json.dumps({key:result[key] for key in ('accepted','gates','qacc_max','COM_acceleration','ground_N','loop_gaps_m','wall_s')}, indent=2))


if __name__ == '__main__':
    main()
