"""Read-only strict compatibility check against the existing accelerated W2L kernel."""
from pathlib import Path
import sys,json
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031'))
import native_kernel as K
meta=json.loads((HERE/'build/runs/Biarticular_net/IMPORT_AUDIT.json').read_text())
out=dict(status='not yet attempted',kernel='virtual_walker_20261031/native_kernel.py',changes_to_source_circuit=False)
try:
 c=K.NativeController(HERE/'Biarticular_net_Standalone.asim',meta['muscles'])
 out.update(status='constructor accepted',neurons=len(c.names),edges=len(c.edges))
except Exception as e:
 out.update(status='existing W2L accelerated controller rejected source',error_type=type(e).__name__,error=str(e),
  scope='No fallback neural circuit or altered source wiring substituted; source-specific binding/parity work required before coupled MuJoCo claims.')
(HERE/'ACCELERATED_PORT_PRECHECK.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
