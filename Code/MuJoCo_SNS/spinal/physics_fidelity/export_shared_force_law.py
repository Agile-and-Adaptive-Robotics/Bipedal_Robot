"""Export a compact opt-in parameter file from the completed static audit."""
import json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
raw=HERE/'runs/audit_20261030/SHARED_FORCE_SURROGATE.json'
r=json.loads(raw.read_text());families={}
for name,v in r['muscles'].items():
 assert v['joint_fit']['heldout_RMSE_N']<v['baseline']['heldout_RMSE_N'],name
 assert v['joint_fit']['maximum_bilateral_force_difference_N']<.01,name
 assert v['original_parameters']['r']['dynamics']==v['original_parameters']['l']['dynamics'],name
 families[name]=dict(actuator_names=[name+'_r',name+'_l'],gainprm=v['joint_gainprm'],biasprm=v['joint_biasprm'],dynprm=v['joint_dynprm'],lengthrange=v['shared_lengthrange'],
  validation=dict(original=v['baseline'],canonical_right=v['canonical_right'],shared=v['joint_fit']))
out=dict(format='aarl_shared_surrogate_force_law_v1',model_family='gait2392_simbody_cvt3_repaired',
 scope='OPT-IN diagnostic: six bilateral muscle families, common static MuJoCo force-law parameters jointly fitted to fresh OpenSim equilibrium forces. Not full-ROM/dynamic/compliant-tendon validation. Other muscles and all geometry/inertia/contact parameters remain outside this file.',
 activations=list(r['activations']),training_poses=[p for p in r['poses'] if not p['heldout']],heldout_poses=[p for p in r['poses'] if p['heldout']],
 source_opensim=r['source_opensim'],source_audit_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),
 loader_contract='For each actuator_names entry, copy gainprm,biasprm,dynprm,lengthrange exactly. Preserve actuator gain/bias/dynamics types, controls, activation bounds, tendon routes, timestep, and integrator. Recompute any cached Fmax or muscle normalization AFTER loading. Reject unknown family/name/shape, nonfinite parameters, or model mismatch. Record file hash in trial provenance. Do not infer a production default.',
 muscles=families)
p=HERE/'shared_force_law_functional_v1.json';p.write_text(json.dumps(out,indent=2));print(p)
