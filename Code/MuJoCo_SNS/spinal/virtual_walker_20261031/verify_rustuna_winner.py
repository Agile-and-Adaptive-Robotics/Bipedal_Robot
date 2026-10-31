"""Fresh-process repeat of the bounded calibration winner."""
import json,subprocess,sys
import numpy as np
from import_virtual_walker import HERE
from assess_trial import assess
best=json.loads((HERE/'runs/rustuna108_ground_home_v1/BEST.json').read_text())['trial'];cfg=json.loads((HERE/'CANONICAL108_NEARFLAT_ANKLE_PHASE.json').read_text());cfg.update(best['params']);path=HERE/'RUSTUNA108_GROUND_WINNER.json';path.write_text(json.dumps(cfg,indent=2));tag='canonical108_ground_rustuna_winner_fresh_v1'
cfg['source']=str(HERE.parents[3]/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer_Standalone.asim');cfg['expected_source_sha256']=best['diagnostics']['source_sha256'];path.write_text(json.dumps(cfg,indent=2))
with (HERE/'runs/rustuna108_winner_fresh.log').open('w') as log:subprocess.run([sys.executable,str(HERE/'run_virtual_walker.py'),'--newer','--config',str(path),'--duration','10','--tag',tag],stdout=log,stderr=subprocess.STDOUT,check=True)
folder=HERE/'runs/W2L_VirtualWalker_Newer'/tag;a=assess(folder);original=np.load(best['diagnostics']['folder']+'/trajectory.npz');fresh=np.load(folder/'trajectory.npz');parity={k:original[k].shape==fresh[k].shape and original[k].dtype==fresh[k].dtype and original[k].tobytes()==fresh[k].tobytes() for k in original.files if original[k].dtype.kind in 'fiu'};assert all(parity.values()),parity
out=dict(config=str(path),selected_trial=best['number'],parameters=best['params'],all_numeric_arrays_byte_exact=True,arrays=parity,assessment=a,result=json.loads((folder/'RESULT.json').read_text()));(folder/'WINNER_FRESH_REPEAT.json').write_text(json.dumps(out,indent=2));(HERE/'runs/RUSTUNA108_WINNER_VERIFICATION.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
