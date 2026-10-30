"""New shared-parameter extension; preserve v1/v2 and disclose route residuals."""
import json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3];RAW=HERE/'runs/audit_20261030'
v1=HERE/'shared_force_law_functional_v1.json';v2=HERE/'shared_force_law_functional_v2.json'
originals={p:p.read_bytes() for p in (v1,v2)};out=json.loads(originals[v2]);accepted=[];rejected={};validation={}
profiles=[('hip_extension_v2','SHARED_FORCE_SURROGATE_HIP_V2.json',('iliacus','lat_gas')),('biarticular_extension_v3','SHARED_FORCE_SURROGATE_BIARTICULAR_V3.json',None)]
for profile,filename,names in profiles:
 path=RAW/filename;r=json.loads(path.read_text());validation[profile]=dict(source_audit_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),activations=r['activations'],training_poses=[p for p in r['poses'] if not p['heldout']],heldout_poses=[p for p in r['poses'] if p['heldout']])
 for name,v in r['muscles'].items():
  if names is not None and name not in names:continue
  errors={}
  for key in ('baseline_tension_N','canonical_right_tension_N','joint_tension_N'):
   values=[];held=[]
   for a,b in zip(v['points'][0::2],v['points'][1::2]):
    assert a['side']=='r' and b['side']=='l' and a['pose']==b['pose'] and a['activation']==b['activation']
    error=(a[key]-b[key])-(a['opensim_tension_N']-b['opensim_tension_N']);values.append(abs(error))
    if a['pose']['heldout']:held.append(abs(error))
   errors[key]=dict(max_source_bilateral_difference_error_N=max(values),heldout_max_source_bilateral_difference_error_N=max(held))
  reasons=[];joint=errors['joint_tension_N']['max_source_bilateral_difference_error_N'];baseline=errors['baseline_tension_N']['max_source_bilateral_difference_error_N']
  if not v['joint_fit']['heldout_RMSE_N']<v['baseline']['heldout_RMSE_N']:reasons.append('Heldout source-force RMSE does not improve')
  if joint>baseline+1e-6:reasons.append('Bilateral force-difference error versus source worsens')
  if joint/v['source_Fmax_N']>.01:reasons.append('Residual force-difference error exceeds1% sourceFmax diagnostic bound')
  if v['original_parameters']['r']['dynamics']!=v['original_parameters']['l']['dynamics']:reasons.append('Activation dynamics differ')
  if reasons:rejected[name]=dict(reasons=reasons,force_difference_validation=errors,fit=v['joint_fit']);continue
  assert name not in out['muscles']
  out['muscles'][name]=dict(actuator_names=[name+'_r',name+'_l'],gainprm=v['joint_gainprm'],biasprm=v['joint_biasprm'],dynprm=v['joint_dynprm'],lengthrange=v['shared_lengthrange'],validation_profile=profile,
   validation=dict(original=v['baseline'],canonical_right=v['canonical_right'],shared=v['joint_fit'],source_bilateral_difference=errors))
  accepted.append(name)
out['variant']='functional_v3';out['base_parameter_file']=v2.name;out['base_parameter_sha256']=hashlib.sha256(originals[v2]).hexdigest()
out['scope']='OPT-IN static shared-PARAMETER diagnostic. Preserve every v2 law, append accepted hip/knee biarticular laws and previously withheld iliacus/lateral-gastrocnemius when source-force error improves. Geometry can leave small delivered-force differences; bilateral difference is compared with source OpenSim, not forced to zero. Functional static data only; no full-ROM/dynamic/compliant-tendon validation.'
out['extension_validation_v3']=dict(profiles=validation,accepted_families=accepted,rejected_families=rejected,
 selection='Improve heldout source-force RMSE and do not worsen maximum bilateral FORCE-DIFFERENCE error versus source; residual <=1% sourceFmax is an explicit diagnostic bound, not a physiological standard.')
p=HERE/'shared_force_law_functional_v3.json';p.write_text(json.dumps(out,indent=2))
assert all(q.read_bytes()==content for q,content in originals.items())
print(p);print('Accepted',accepted);print('Rejected',rejected)
