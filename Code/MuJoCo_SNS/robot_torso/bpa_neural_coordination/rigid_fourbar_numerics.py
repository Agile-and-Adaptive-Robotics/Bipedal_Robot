"""Explicit opt-in engineering numericalapproximation, not physicalcalibration."""
import json,hashlib
from pathlib import Path
import numpy as np
import mujoco
DEFAULT_RECORD=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\cad_bpa\pair_neural_recruitment\normalized_unclamped\fourbar_numerical_fidelity')

def apply(model,enabled=False,record=DEFAULT_RECORD,iterations=500,timestep=None):
 if not enabled:return dict(enabled=False,changed=False)
 record=Path(record);forward=json.loads((record/'FIXTURE_FORWARD_AUDIT.json').read_text());signed=json.loads((record/'SIGNED_BILATERAL_CHECK.json').read_text());normal=json.loads((record/'NORMAL_LOAD_FORWARD_AUDIT.json').read_text());assert forward['selected_fixture_gate'] and forward['time_convergence_pass'] and signed['gate'] and normal['parent_05mm_gate'] and normal['half_quarter_converged'],'Numericalfixturevalidation is not passed';selection=json.loads((record/'DERIVED_VARIANT.json').read_text());protected={name:getattr(model,name).copy() for name in ('body_mass','body_inertia','body_ipos','body_iquat','dof_damping','dof_armature','jnt_range')};before={};ids=[]
 for side in ('r','l'):
  idx=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_EQUALITY,'knee_loop_'+side);assert idx>=0;ids.append(idx);before[side]=dict(solref=model.eq_solref[idx].tolist(),solimp=model.eq_solimp[idx].tolist());model.eq_solref[idx]=[selection['timeconst_s'],1.];model.eq_solimp[idx,:2]=selection['constant_impedance']
 old_dt=float(model.opt.timestep);old_iterations=int(model.opt.iterations);model.opt.iterations=iterations
 if timestep is not None:model.opt.timestep=timestep
 assert model.opt.timestep<=selection['timeconst_s']/2.+1e-15,'refsafe numericaltimeconstant condition violated'
 for name,array in protected.items():assert np.array_equal(array,getattr(model,name)),name+'changed'
 return dict(enabled=True,scope='Engineeringnumericalrigid-linkapproximation, sourcegeometry/mass/damping unchanged; no hardwarestiffness or robotstandingvalidation',old=before,new={side:dict(solref=model.eq_solref[idx].tolist(),solimp=model.eq_solimp[idx].tolist()) for side,idx in zip(('r','l'),ids)},old_dt=old_dt,new_dt=float(model.opt.timestep),old_iterations=old_iterations,new_iterations=int(model.opt.iterations),protected_arrays_exact=True,fixture_record=str(record),selection_sha256=hashlib.sha256((record/'DERIVED_VARIANT.json').read_bytes()).hexdigest())
