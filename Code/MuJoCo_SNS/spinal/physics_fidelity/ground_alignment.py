"""Opt-in foot-plane pose alignment before support-rig references are built."""
import numpy as np
import mujoco
from scipy.optimize import brentq
from .contact_standing_solver import _support_regions


def align_ground_pose(model,data,*,flatten_ankles=True,target_gap_m=-.001,ankle_bounds_deg=(-20.,20.)):
    """Mutate pose only; return exact alignment/reference metadata.

    Model MUST already have the actual repaired joint axes. A raw converted
    model with pre-repair hip axes is invalid; use align_via_repaired_probe.
    Call AFTER selecting hip/knee/etc pose and BEFORE capture_pose/apply_harness.
    A previously constructed vertical rig must be rebased by the returned
    vertical-coordinate shift; leaving its old spring reference is invalid.
    Target gap is explicit signed mesh-plane distance (negative penetration),
    not a calibrated contact force or a guarantee of static equilibrium.
    """
    from mechanics import CoupledArms
    ar=CoupledArms(model);ar.project(data);mujoco.mj_forward(model,data)
    before,normal=_support_regions(model,data);ankles=[model.joint('ankle_angle_'+s).qposadr[0] for s in ('r','l')]
    old_angles=data.qpos[ankles].copy();chosen=float(old_angles.mean())
    if flatten_ankles:
        def mismatch(angle):
            data.qpos[ankles]=angle;ar.project(data);mujoco.mj_forward(model,data)
            points,_=_support_regions(model,data)
            return .5*((points[1]['gap_m']-points[0]['gap_m'])+(points[3]['gap_m']-points[2]['gap_m']))
        lo,hi=np.radians(ankle_bounds_deg);flo,fhi=mismatch(lo),mismatch(hi)
        if flo*fhi>0.:
            data.qpos[ankles]=old_angles;ar.project(data);mujoco.mj_forward(model,data)
            raise ValueError('No common ankle flattening root within declared bounds')
        chosen=float(brentq(mismatch,lo,hi,xtol=1e-11));mismatch(chosen)
    points,normal=_support_regions(model,data);gap=min(p['gap_m'] for p in points)
    candidates=[model.joint('pelvis_'+s).id for s in ('tx','ty','tz')]
    j=max(candidates,key=lambda i:abs(float(data.xaxis[i]@normal)));projection=float(data.xaxis[j]@normal)
    if abs(projection)<.8:raise ValueError('No pelvis slide aligned with ground normal')
    qadr=int(model.jnt_qposadr[j]);shift=(target_gap_m-gap)/projection;data.qpos[qadr]+=shift
    ar.project(data);mujoco.mj_forward(model,data);after,_=_support_regions(model,data)
    return dict(scope='Explicit initial pose alignment; not a hold/balance/walking result.',flatten_ankles=flatten_ankles,
        target_gap_m=target_gap_m,initial_region_gaps_m=[r['gap_m'] for r in before],final_region_gaps_m=[r['gap_m'] for r in after],
        original_native_ankle_deg=np.degrees(old_angles).tolist(),common_native_ankle_deg=float(np.degrees(chosen)),
        vertical_joint=model.joint(j).name,vertical_qpos_address=qadr,vertical_coordinate_shift_m=float(shift),
        world_normal_displacement_m=float(shift*projection),
        rig_reference_contract='Construct/capture the support rig AFTER this pose. If aligning an already-built rig, shift its vertical qpos_spring by exactly vertical_coordinate_shift_m too; do not silently retain the old reference.')


def align_via_repaired_probe(model,data,harness_builder,**alignment_options):
    """Align a native-coordinate pose held in the raw pre-repair model.

    harness_builder must be the captured ORIGINAL runner.apply_harness, not a
    recursively wrapped pose hook. Final actual rig must be built afterward.
    """
    import runner as R
    pose=R.capture_pose(model,data)
    probe=harness_builder(model,data,kxy=1500.,no_ground=False)
    pd=mujoco.MjData(probe);R.seed_pose(probe,pd,pose);mujoco.mj_forward(probe,pd)
    result=align_ground_pose(probe,pd,**alignment_options)
    R.seed_pose(model,data,R.capture_pose(probe,pd));mujoco.mj_forward(model,data)
    result['pipeline']='Native q values -> repaired probe -> aligned native q values copied back -> capture -> final repaired rig/ref construction.'
    result['probe_raw_hip_axis']=[model.joint('hip_flexion_r').axis.tolist(),probe.joint('hip_flexion_r').axis.tolist()]
    return result
