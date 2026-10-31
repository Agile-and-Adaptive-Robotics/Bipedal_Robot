"""Fresh CLI repeat of lower-hip candidate and unchanged default path."""
import hashlib,json,subprocess,sys
import numpy as np
from import_virtual_walker import HERE
from assess_trial import assess
original=HERE/'RUSTUNA108_GROUND_WINNER.json';original_bytes=original.read_bytes();cfg=json.loads(original_bytes);cfg['hip_pf_mn_scale']=.6;config=HERE/'GENTLE_HIP_0P6.json';config.write_text(json.dumps(cfg,indent=2));base=HERE/'runs/W2L_VirtualWalker_Newer';rows=[]
for label,path,reference in [('gentle_hip_0p6_fresh',config,'canonical108_ground_gentle_hip_0p6_v1'),('hip_default_parity',original,'canonical108_ground_rustuna_winner_fresh_v1')]:
 tag='canonical108_ground_'+label+'_v1'
 with (HERE/'runs'/(label+'.log')).open('w') as log:subprocess.run([sys.executable,str(HERE/'run_virtual_walker.py'),'--newer','--config',str(path),'--duration','10','--tag',tag],stdout=log,stderr=subprocess.STDOUT,check=True)
 folder=base/tag;a=assess(folder);old=np.load(base/reference/'trajectory.npz');new=np.load(folder/'trajectory.npz');parity={k:old[k].shape==new[k].shape and old[k].dtype==new[k].dtype and old[k].tobytes()==new[k].tobytes() for k in old.files if old[k].dtype.kind in 'fiu'};assert all(parity.values()),parity;r=json.loads((folder/'RESULT.json').read_text());out=dict(tag=tag,config=str(path),config_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),all_numeric_arrays_byte_identical=True,arrays=parity,result=r,assessment=a);(folder/'FRESH_HIP_VERIFICATION.json').write_text(json.dumps(out,indent=2));rows.append(out)
assert original.read_bytes()==original_bytes;(HERE/'runs/GENTLE_HIP_FRESH_VERIFICATION.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(original_winner_unchanged=True,fresh_gentle_byte_exact=True,default1_byte_exact=True),indent=2))
