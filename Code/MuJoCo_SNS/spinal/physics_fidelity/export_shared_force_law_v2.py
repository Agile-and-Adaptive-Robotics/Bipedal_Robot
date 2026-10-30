"""Append accepted hip/gastrocnemius static laws without changing v1."""
import json,hashlib,copy
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
base=HERE/'shared_force_law_functional_v1.json';original=base.read_bytes();out=json.loads(original)
raw=HERE/'runs/audit_20261030/SHARED_FORCE_SURROGATE_HIP_V2.json';r=json.loads(raw.read_text());accepted=[];rejected={}
for name,v in r['muscles'].items():
 reasons=[]
 if not v['joint_fit']['heldout_RMSE_N']<v['baseline']['heldout_RMSE_N']:reasons.append('Heldout error does not improve')
 if v['joint_fit']['maximum_bilateral_force_difference_N']>=.01:reasons.append('Bilateral force residual exceeds0.01N')
 if v['original_parameters']['r']['dynamics']!=v['original_parameters']['l']['dynamics']:reasons.append('Activation dynamics differ')
 if reasons:rejected[name]=reasons;continue
 assert name not in out['muscles']
 out['muscles'][name]=dict(actuator_names=[name+'_r',name+'_l'],gainprm=v['joint_gainprm'],biasprm=v['joint_biasprm'],dynprm=v['joint_dynprm'],lengthrange=v['shared_lengthrange'],validation_profile='hip_extension_v2',
  validation=dict(original=v['baseline'],canonical_right=v['canonical_right'],shared=v['joint_fit']))
 accepted.append(name)
out['variant']='functional_v2'
out['scope']='OPT-IN static shared-force diagnostic: six v1 families preserved exactly, plus accepted hip/gastrocnemius extensions fitted on expanded functional poses. Not full-ROM/dynamic/compliant-tendon validation. Other muscle laws and all geometry/inertia/contact parameters remain outside this file.'
out['base_parameter_file']=base.name;out['base_parameter_sha256']=hashlib.sha256(original).hexdigest()
out['extension_validation']=dict(source_audit_sha256=hashlib.sha256(raw.read_bytes()).hexdigest(),training_poses=[p for p in r['poses'] if not p['heldout']],heldout_poses=[p for p in r['poses'] if p['heldout']],accepted_families=accepted,rejected_families=rejected)
p=HERE/'shared_force_law_functional_v2.json';p.write_text(json.dumps(out,indent=2))
assert base.read_bytes()==original
print(p);print('Accepted',accepted,'Rejected',rejected)
