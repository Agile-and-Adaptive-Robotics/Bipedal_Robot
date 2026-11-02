"""One fixed, explicitly distinct5% coactivation fixture; not a fitted winner."""
import sys,json,argparse,hashlib
from pathlib import Path
import numpy as np
import cad_mixed_pressure_experiment as C
import cad_mixed_pressure_live as L
OUT=C.OUT.parent/'mixed_pressure_coactivation05';FLOOR=.05
C.OUT=OUT;L.OUT=OUT
def register():
 OUT.mkdir(exist_ok=True);p=OUT/'COACTIVATION_REGISTRATION.json';d=dict(machine='easteregg2',activation_floor=FLOOR,scope='Existing original hybrid controller activation-floor option, fixed5% across54routes;2pressure-BPAs+52referenceforces in mixed comparison. Not autonomous SNScoordination, all-BPAstanding, or hardwarecalibration.',reason='Original0floor mixed case completes but delivers zero actual tension through both pressure routes; it remains a negative actuation result.',source_pose=1069,nominal_harness_bw=.9,duration_s=2.,same_BPA_geometry=True,no_search=True,stop_on_domain_failure_or_hyperextension=True,zero_floor_results_preserved=True,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(C.__file__),Path(L.__file__))})
 if p.exists():assert json.loads(p.read_text())==d,'Frozen5%fixture changed'
 else:p.write_text(json.dumps(d,indent=2))
 return d
def audit():
 register();original=C.R.allocate
 def fixed(*args,**kw):kw['activation_floor']=FLOOR;return original(*args,**kw)
 C.R.allocate=fixed
 try:return C.audit()
 finally:C.R.allocate=original
def live(mode):
 register();assert json.loads((OUT/'STATIC_AUDIT.json').read_text())['ready_for_live'];original=L.compile_opt_in
 def with_floor():
  fn,parser,code,edits=original();parser.set_defaults(activation_floor=FLOOR);return fn,parser,code,edits
 L.compile_opt_in=with_floor
 return L.run(mode)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['reference','mixed']);a=p.parse_args()
 if a.mode:live(a.mode)
 else:audit()
