"""Static, equality-reduced mechanical descriptors for Sutton-style scaling.

No time integration, optimization, harness tuning, or controller fitting.
The pelvis is held at each reference pose for this internal-coordinate map.
"""
import argparse
import json
import os
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import mujoco
import numpy as np
from syn6_torque_span import MODEL, REFERENCE, SPINAL, REPO, CoupledArms, JOINTS, sha

OSIM = REPO / 'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim'
FORCE = SPINAL / 'physics_fidelity/shared_force_law_functional_v3.json'


def load_variant(variant):
    m = mujoco.MjModel.from_xml_path(str(MODEL))
    if variant == 'raw_cvt3':
        return m, {'patch': None, 'accepted_standing_configuration_matched': False}
    if variant != 'runtime_patch_v3':
        raise ValueError('Unknown explicit physics variant: ' + variant)
    sys.path[:0] = [str(SPINAL), str(SPINAL / 'w2l_cpg')]
    import runner as R
    d = mujoco.MjData(m)
    R.apply_start_pose(m, d)
    text = R.patch_xml(MODEL.read_text(encoding='utf-8'), m, d)
    directory = MODEL.parent.resolve()
    fd, temporary = tempfile.mkstemp(prefix='_fsa_mechanical_audit_', suffix='.xml', dir=directory)
    target = Path(temporary).resolve()
    if target.parent != directory or not target.name.startswith('_fsa_mechanical_audit_'):
        os.close(fd)
        raise ValueError('Temporary model path escaped intended directory')
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as handle:
            handle.write(text)
        m = mujoco.MjModel.from_xml_path(str(target))
    finally:
        target.unlink(missing_ok=True)
    force = json.loads(FORCE.read_text())
    for row in force['muscles'].values():
        for name in row['actuator_names']:
            aid = m.actuator(name).id
            for key in ('gainprm', 'biasprm', 'dynprm', 'lengthrange'):
                getattr(m, 'actuator_' + key)[aid] = row[key]
    patched_tags = {j.attrib['name']:j.attrib for j in ET.fromstring(text).findall('.//worldbody//joint') if 'name' in j.attrib}
    return m, {'patch': 'runner.patch_xml with all signature defaults, followed by F.run v3 muscle parameter copy',
               'actual_patched_joint_damping_attributes':{name:patched_tags[name].get('damping') for name in JOINTS},
               'damping_realization_note':'Only existing damping attributes are replaced in patch_xml; absent explicit attributes remain inherited. Use compiled DOF damping, not intended tuple values.',
               'shared_force_muscle_families': sorted(force['muscles']),
               'accepted_standing_configuration_matched': False,
               'not_replicated_standing_overrides': ['Corrected saved qfull start pose', 'Knee zero-extension margin/stop override',
                  'F.run vertical rig stiffness times0.01 and damping times0.42', 'Any supplied per-axis rig scale or other run-specific override'],
               'spring_reference_source': 'Raw model after runner.apply_start_pose, before compiling patch; constant for all four static samples'}


def body_ledger(m, arms):
    source = ET.parse(OSIM).getroot().find('.//BodySet/objects')
    bodies = {b.attrib['name']: float(b.findtext('mass')) for b in source if b.tag == 'Body'}
    dependent = set(arms.dependent.tolist())
    rows = []
    for bid in range(1, m.nbody):
        name = mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_BODY, bid)
        joints = [j for j in range(m.njnt) if m.jnt_bodyid[j] == bid]
        follower = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_JOINT, j) for j in joints if m.jnt_dofadr[j] in dependent]
        rows.append(dict(name=name, model_mass_kg=float(m.body_mass[bid]), source_anatomical_body=name in bodies,
                         fixed_to_world=bool(m.body_weldid[bid]==0),
                         source_BodySet_mass_kg=bodies.get(name), equality_dependent_joints=follower))
    return dict(source_osim=str(OSIM), source_BodySet_names=sorted(bodies), source_BodySet_total_mass_kg=sum(bodies.values()),
                moving_body_mass_kg=sum(r['model_mass_kg'] for r in rows if not r['fixed_to_world']),
                fixed_to_world_mass_kg=sum(r['model_mass_kg'] for r in rows if r['fixed_to_world']),
                matched_anatomical_mass_kg=sum(r['model_mass_kg'] for r in rows if r['source_anatomical_body']),
                non_BodySet_auxiliary_or_unmatched_mass_kg=sum(r['model_mass_kg'] for r in rows if not r['source_anatomical_body']),
                dependent_joint_anatomical_mass_kg=sum(r['model_mass_kg'] for r in rows if r['source_anatomical_body'] and r['equality_dependent_joints']),
                dependent_joint_auxiliary_or_unmatched_mass_kg=sum(r['model_mass_kg'] for r in rows if not r['source_anatomical_body'] and r['equality_dependent_joints']),
                interpretation='A dependent joint is not evidence of an artificial body: tibia_r/l are genuine source BodySet segments.', bodies=rows)


def audit(variant='raw_cvt3'):
    inputs = (MODEL, REFERENCE, OSIM, FORCE, SPINAL/'runner.py', SPINAL/'functional_gait2392/mechanics.py')
    frozen = {str(p): sha(p) for p in inputs}
    m, construction = load_variant(variant)
    d = mujoco.MjData(m)
    arms = CoupledArms(m)
    with np.load(REFERENCE, allow_pickle=False) as z:
        names = [str(x) for x in z["coordinate_names"]]
        values, phase = z["r_mean"].copy(), z["phase"].copy()
    all_jids = [mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, name) for name in names]
    if min(all_jids) < 0 or m.nq != m.nv:
        raise ValueError("This explicit adapter requires scalar anatomical coordinates")
    full_mass = np.zeros((m.nv, m.nv))

    def terms(qpos, velocity=None):
        d.qpos[:] = qpos
        d.qvel[:] = 0.
        d.ctrl[:] = 0.
        d.act[:] = 0.
        arms.project(d)
        if velocity is not None:
            arms.update(d)
            d.qvel[:] = arms.E @ velocity
        mujoco.mj_forward(m, d)
        arms.update(d)
        E = arms.E.copy()
        gravity = E.T @ d.qfrc_bias
        # Native muscle bias force at zero activation and coordinate/tendon
        # springs are included; contacts/constraint reactions are not included.
        restoring = gravity - E.T @ (d.qfrc_passive + d.qfrc_actuator)
        return E, gravity.copy(), restoring.copy()

    samples = []
    for percent in (0., 25., 50., 75.):
        index = int(np.argmin(np.abs(phase - percent)))
        pose = m.qpos0.copy()
        for jid, name, value in zip(all_jids, names, values[index]):
            pose[m.jnt_qposadr[jid]] = value if name in ("pelvis_tx", "pelvis_ty", "pelvis_tz") else np.radians(value)
        for name in ('knee_angle_r', 'knee_angle_l'):
            if pose[m.joint(name).qposadr[0]] > 0:
                raise ValueError('Reference knee hyperextension: do not silently clip')
        E, gravity, restoring = terms(pose)
        pose = d.qpos.copy()
        mujoco.mj_fullM(m, full_mass, d.qM)
        reduced_mass = E.T @ full_mass @ E
        coordinate_damping = E.T @ (m.dof_damping[:, None] * E)
        damping = np.zeros((17, 17))
        dv = np.zeros(17)
        for column in range(17):
            dv[:] = 0.; dv[column] = 1e-5
            _, _, rminus = terms(pose, -dv)
            _, _, rplus = terms(pose, dv)
            damping[:, column] = (rplus-rminus)/2e-5
        armature = E.T @ (m.dof_armature[:, None] * E)
        following_armature = E[arms.dependent].T @ (m.dof_armature[arms.dependent, None] * E[arms.dependent])
        eig = np.linalg.eigvalsh(reduced_mass)
        if eig[0] <= 0:
            raise ValueError("Reduced inertia is not positive definite")
        Kgravity, Ktotal = np.zeros((17, 17)), np.zeros((17, 17))
        increment = 1e-5
        for column, address in enumerate(arms.addresses):
            minus = pose.copy(); minus[address] -= increment
            plus = pose.copy(); plus[address] += increment
            _, g0, r0 = terms(minus)
            _, g1, r1 = terms(plus)
            Kgravity[:, column] = (g1 - g0) / (2 * increment)
            Ktotal[:, column] = (r1 - r0) / (2 * increment)
        # Internal rotations all use radians/Nm. This is the local torque
        # impedance for the stated fixed-pelvis, zero-activation condition.
        frequencies = {}
        for label, frequency in (("walking_reference", 1 / 1.2333), ("fast_cached_stride", 1 / .559), ("fast_ankle_angle_band", 4 / .559)):
            omega = 2 * np.pi * frequency
            Z = -omega**2 * reduced_mass + 1j * omega * damping + Ktotal
            frequencies[label] = {"frequency_hz": frequency,
                "one_degree_single_coordinate_torque_norm_Nm": np.linalg.norm(Z * np.radians(1.), axis=0).tolist()}
        samples.append({"phase_percent": float(phase[index]), "mass_kg_m2": reduced_mass.tolist(),
            "joint_damping_Nm_s_rad": damping.tolist(), "gravity_stiffness_Nm_rad": Kgravity.tolist(),
            "projected_coordinate_viscous_damping_Nm_s_rad": coordinate_damping.tolist(),
            "native_velocity_derivative_minus_coordinate_damping_max": float(np.max(np.abs(damping-coordinate_damping))),
            "total_zero_activation_restoring_stiffness_Nm_rad": Ktotal.tolist(),
            "mass_eigenvalues": eig.tolist(), "diagonal_inertia_other_internal_coordinates_held_kg_m2": np.diag(reduced_mass).tolist(),
            "inertia_other_internal_coordinates_free_kg_m2": (1 / np.diag(np.linalg.inv(reduced_mass))).tolist(),
            "armature_diagonal_kg_m2": np.diag(armature).tolist(),
            "follower_armature_diagonal_kg_m2": np.diag(following_armature).tolist(),
            "gravity_Hessian_asymmetry_Nm_rad": float(np.max(np.abs(Kgravity - Kgravity.T))),
            "restoring_Hessian_asymmetry_Nm_rad": float(np.max(np.abs(Ktotal - Ktotal.T))),
            "gravity_torque_Nm": gravity.tolist(), "zero_activation_restoring_torque_Nm": restoring.tolist(),
            "source_reference_index": index,
            "internal_angles_deg": dict(zip(JOINTS,np.degrees(pose[arms.addresses]).tolist())),
            "fixed_pelvis_positions_native": {n:float(pose[m.joint(n).qposadr[0]]) for n in names if n.startswith('pelvis_')},
            "frequency_context": frequencies})
    follower_joints = [i for i in range(m.njnt) if m.jnt_dofadr[i] in set(arms.dependent)]
    follower_bodies = sorted(set(int(m.jnt_bodyid[i]) for i in follower_joints))
    if frozen != {str(p): sha(p) for p in inputs}:
        raise RuntimeError("Input changed during static mechanical audit")
    return {"machine": "easteregg2", "variant":variant,"construction":construction,"source_sha256": frozen, "no_time_integration_or_training": True,
        "body_mass_ledger":body_ledger(m,arms),
        "software_nq":m.nq,"software_nv":m.nv,"floating_base_anatomical_coordinates":6,"internal_anatomical_coordinates":17,
        "reference_scope":"Same four0/25/50/75 percent samples of one cached old walking cycle; no running or equilibrium validation",
        "coordinates": list(JOINTS), "compiled_model_total_mass_kg": float(m.body_mass.sum()),
        "follower_bodies": len(follower_bodies), "follower_body_mass_kg": float(m.body_mass[follower_bodies].sum()),
        "reduction": "Fixed pelvis; equality-aware17internal coordinates. Mred=E.T*M*E; Dred=d[E.T*(qbias-qpassive-qactuator)]/dv at v=0; K=d[E.T*(qbias-qpassive-qactuator)]/dq at activation0",
        "joint_spring_damping_ledger":[dict(name=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,int(j)),stiffness_Nm_rad=float(m.jnt_stiffness[j]),damping_Nm_s_rad=float(m.dof_damping[m.jnt_dofadr[j]]),spring_reference_rad=float(m.qpos_spring[m.jnt_qposadr[j]])) for j in arms.jids],
        "scaling_use": "Compare actual inertial/damping/restoring torque terms at a declared frequency; do not multiply neural clocks by mass alone",
        "limitations": ["Native converted source force law, not source-OpenSim or compliant-tendon/BPA equivalence validation",
            "No harness/contact stiffness, equilibrium solve, actuator activation dynamics, delay, or active-feedback stiffness included",
            "Gravity/restoring stiffness can be negative; no positive scalar natural frequency is inferred from it",
            "Anatomical tibiae have dependent knee-translation joints; their real source mass is separated from non-BodySet auxiliary bodies",
            "Constraint/contact reactions and their effective stiffness/damping are excluded; explicit internal rig springs ARE included",
            "Raw conversion is not current runtime physics; default patched-v3 is not asserted to reproduce saved standing overrides",
            "This internal subsystem descriptor is not a closed-loop standing or walking stability calculation"],
        "samples": samples}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--variant',choices=('raw_cvt3','runtime_patch_v3','both'),default='raw_cvt3')
    args = parser.parse_args()
    if args.variant=='both' and not args.output_dir:parser.error('--variant both requires --output-dir')
    if not args.output_dir and not args.output:parser.error('Supply --output or --output-dir')
    results=[]
    for variant in (('raw_cvt3','runtime_patch_v3') if args.variant=='both' else (args.variant,)):
        result = audit(variant);results.append(result)
        output=args.output_dir/('GAIT_MECHANICAL_SCALING_'+variant+'.json') if args.output_dir else args.output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(result, indent=2), encoding='utf-8')
        print(json.dumps({"report": str(output), "mass_kg": result["compiled_model_total_mass_kg"],
        "follower_mass_kg": result["follower_body_mass_kg"],
        "minimum_inertia_eigenvalue": min(row["mass_eigenvalues"][0] for row in result["samples"]),
        "max_gravity_hessian_asymmetry": max(row["gravity_Hessian_asymmetry_Nm_rad"] for row in result["samples"])}))
    if len(results)==2:
        lines=['# Gait2392 static mechanical scaling comparison','','Fixed pelvis,17 equality-reduced internal coordinates; no integration or fitting.','',
               '| Variant | Compiled mass sum kg | Moving body mass kg | Matched anatomical kg | Other body mass kg | Knee R/L held-coordinate inertia kg m² at0% | Follower armature knee R/L kg m² |','|---|---:|---:|---:|---:|---|---|']
        for r in results:
            b=r['body_mass_ledger'];s=r['samples'][0];ix=[JOINTS.index('knee_angle_'+side) for side in ('r','l')]
            lines.append(f"| {r['variant']} | {r['compiled_model_total_mass_kg']:.9f} | {b['moving_body_mass_kg']:.9f} | {b['matched_anatomical_mass_kg']:.9f} | {b['non_BodySet_auxiliary_or_unmatched_mass_kg']:.9f} | "+'/'.join(f"{s['diagonal_inertia_other_internal_coordinates_held_kg_m2'][i]:.9g}" for i in ix)+' | '+'/'.join(f"{s['follower_armature_diagonal_kg_m2'][i]:.9g}" for i in ix)+' |')
        lines+=['','Tibia_r/l are genuine OpenSim BodySet segments of 3.7075 kg each, despite carrying dependent knee translations. They are not auxiliary pathpoint mass.',
                '', 'Runtime boundmass assigns 0.1 kg to each of 32 moving pathpoint bodies (3.2 kg total), and 0.1 kg to fixed ground. The fixed ground mass is included in the compiled body-mass sum but contributes no moving inertia.',
                '', 'Compiled damping mismatch: patch_xml only replaces existing damping attributes. These internal joint tags inherit damping, so lumbar/hip-rotation damping remains 0.05 rather than the intended 25 Nm s/rad, and ankle/subtalar/MTP/hip-adduction damping remains 0.05 rather than the intended 3. This audit records the actual compiled damping, makes no repair, and preserves the source hash. The runtime JSON includes actual patched attribute presence and compiled joint values.',
                '', 'The patched variant includes runner.patch_xml default boundmass, follower armature, axis repairs, coordinate springs/damping, muscle pruning, and v3 shared static force-law parameters. It does not replicate the corrected standing start pose, protective knee override, or reduced vertical rig. Internal lumbar/hip-rotation rig springs are included; floating-pelvis rig motion, native contact/constraint stiffness, activation/fiber dynamics, and feedback are excluded.',
                '', 'K is evaluated from projected native gravity minus coordinate/tendon passive and zero-activation muscle forces, including pose-dependent equality projection. D is independently differentiated with velocity mapped through that projection. Negative stiffness or coupled inertia cannot be replaced by a positive scalar mass/time-scale heuristic. These are local torque descriptors at four old walking postures, not equilibrium or closed-loop stability evidence.']
        (args.output_dir/'GAIT_MECHANICAL_SCALING_COMPARISON.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
