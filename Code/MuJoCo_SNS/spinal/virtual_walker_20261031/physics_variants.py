"""Explicit opt-in runtime mechanical changes; source XML and geometry stay intact."""
import math
import numpy as np
import mujoco

def apply_segment_mass_scale(model,data,scale=1.0):
 """Scale bilateral limb density at fixed shape, including its principal inertia.

 Call after model compilation and before first forward/measurement. The default
 returns without mj_setConst or any array assignment, preserving the baseline.
 """
 scale=float(scale)
 if not math.isfinite(scale) or scale<=0:raise ValueError('segment_mass_scale must be finite and positive')
 names=[segment+'_'+side for side in ('L','R') for segment in ('femur','tibia','foot','toe')]
 ids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,name) for name in names]
 if any(i<0 for i in ids):raise ValueError('Mass variation requires all eight named bilateral limb bodies')
 ledger=[dict(body=name,body_id=int(i),mass_before_kg=float(model.body_mass[i]),inertia_before_kg_m2=model.body_inertia[i].copy().tolist(),mass_after_kg=float(model.body_mass[i]*scale),inertia_after_kg_m2=(model.body_inertia[i]*scale).tolist()) for name,i in zip(names,ids)]
 before=float(np.sum(model.body_mass[ids]))
 if scale!=1.0:
  model.body_mass[ids]*=scale;model.body_inertia[ids]*=scale
  mujoco.mj_setConst(model,data)
 return dict(segment_mass_scale=scale,applied=scale!=1.,body_ledger=ledger,selected_mass_before_kg=before,selected_mass_after_kg=before*scale,
  geometry_changed=False,source_modified=False,scope='Fixed-geometry bilateral limb-density perturbation only; pelvis, apparatus, contact-only bodies, neural parameters, and tendon paths unchanged')
