"""Audit ideal-reference activation semantics without integrating mechanics."""
import json
import numpy as np
import cad_rigid_stance as C

OUT = C.OUT / 'zero_floor_capacity_case'


def audit():
    OUT.mkdir(exist_ok=True)
    b = C.Native()
    source = C.F.B.S.P.C.ROBOT / 'ah001_neural.xml'
    model = C.mujoco.MjModel.from_xml_path(str(source))
    data = C.mujoco.MjData(model)
    state = C.P.load(C.OUT / 'REFINED_PRELOAD_STATE.npz')
    capacities = -model.actuator_gainprm[:, 0]
    maps = []
    for name, qpos in [('exact_closed_geometry', b.base), ('refined_geometry', state['qpos'])]:
        data.qpos[:] = qpos
        data.qvel[:] = 0
        C.mujoco.mj_forward(model, data)
        for activation in (0., .05, .1):
            data.act[:] = activation
            data.ctrl[:] = activation
            C.mujoco.mj_fwdActuation(model, data)
            maps.append(dict(pose=name, activation=activation,
                             tension_N=(-data.actuator_force).tolist(),
                             map_error_N=float(np.max(np.abs(data.actuator_force + capacities * activation)))))
    assert all(row['map_error_N'] < 1e-12 for row in maps)
    floor = .05 * capacities
    floor[b.s.ids] = 0
    b.d.qpos[:] = state['qpos']
    b.d.qvel[:] = 0
    b.s.dynamic.update(b.d)
    generalized = -b.s.dynamic.length_jacobian[b.s.tids].T @ floor
    result = dict(machine='easteregg2', source=str(source), native_actuator_count=model.nu,
                  all_fixed_gain=bool(np.all(model.actuator_gaintype == int(C.mujoco.mjtGain.mjGAIN_FIXED))),
                  all_no_bias=bool(np.all(model.actuator_biastype == int(C.mujoco.mjtBias.mjBIAS_NONE))),
                  all_first_order_filter=bool(np.all(model.actuator_dyntype == int(C.mujoco.mjtDyn.mjDYN_FILTER))),
                  time_constants_s=model.actuator_dynprm[:, 0].tolist(), static_map=maps,
                  source_definition='da/dt=(u-a)/0.03; tensile_reference=F0*a. No force-length, force-velocity, passive-fiber, or tendon law.',
                  physiological_5percent_activation_map_available=False,
                  reason='Only human-named peak capacities were assigned to these custom CAD paths; route-specific muscle-law parameters are absent.',
                  prior_other52_floor_N=floor.tolist(), sum_internal_reference_tensions_N=float(floor.sum()),
                  largest_single_floor_N=float(floor.max()), source_joint_names=[model.joint(int(j)).name for j in model.dof_jntid],
                  floor_generalized_force_at_refined_pose=generalized.tolist(),
                  interpretation='The sum is internal route tension, not an external support load. Five percent of this ideal state equals five percent of peak reference tension, not validated physiological coactivation.',
                  zero_floor_fixture_definition='T_i in [0,F0_i] for all54, no modeled passive tension, no SNS or physiological-coactivation claim. Mathematical/mechanical force feasibility only.',
                  mechanical_integration_performed=False)
    (OUT / 'ACTIVATION_SEMANTICS_AUDIT.json').write_text(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    result = audit()
    print(json.dumps({key: result[key] for key in ('all_fixed_gain', 'all_no_bias', 'all_first_order_filter', 'physiological_5percent_activation_map_available', 'sum_internal_reference_tensions_N', 'largest_single_floor_N')}, indent=2))
