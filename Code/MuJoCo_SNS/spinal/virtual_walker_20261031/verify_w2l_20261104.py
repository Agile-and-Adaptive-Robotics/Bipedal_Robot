"""Fresh-process physical repeat and saved-state preview of the bounded study."""
import argparse,json,contextlib
from pathlib import Path
import numpy as np
from run_virtual_walker import run
from import_virtual_walker import HERE
from assess_trial import assess
from render_trial import render

def main():
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();summary=json.loads((a.root/'SUMMARY.json').read_text());best=summary['best'];cfg=json.loads((HERE/'CANONICAL108_NEARFLAT_ANKLE_PHASE.json').read_text());cfg.update(best['params'])
 with(a.root/'fresh_repeat.log').open('w',encoding='utf-8')as log,contextlib.redirect_stdout(log):
  result=run(newer=True,duration=10.,tag='repeat',output_root=a.root/'fresh_repeat',**cfg)
  folder=a.root/'fresh_repeat/W2L_VirtualWalker_Newer/repeat';assessment=assess(folder)
  old=np.load(Path(best['diagnostics']['folder'])/'trajectory.npz');new=np.load(folder/'trajectory.npz');parity={key:old[key].shape==new[key].shape and old[key].tobytes()==new[key].tobytes()for key in old.files if old[key].dtype.kind in 'fiu'}
  proof=dict(parameters=best['params'],selected_trial=best['number'],all_numeric_arrays_byte_exact=all(parity.values()),arrays=parity,assessment=assessment,result=result,measured_assistance=best['diagnostics']['support_components'],scope='Li-Deng-derived newer W2L108; fixed apparatus and muscle physics; bilateral scalar recruitment and clock calibration; no free-walking claim')
  (a.root/'FRESH_VERIFICATION.json').write_text(json.dumps(proof,indent=2))
  if not all(parity.values()):raise RuntimeError('Fresh trajectory mismatch')
  render(folder,start=5.,stop=8.2,title='Li-Deng-derived W2L108 | Rustuna calibration',detail_caption=f"Upward assistance after 1 s: {100*best['diagnostics']['upward_assistance_BW']:.1f}% BW | recorded replay")
 print(json.dumps(dict(repeated=all(parity.values()),passed=assessment['alternating_supported_stepping_gate'],folder=str(folder))))
if __name__=='__main__':main()
