"""Exercise the real evaluator and runner without stepping physics.

Every conductance dimension must alter the compiled matrices. Sensor
encoding parameters are separately reported as runtime (not synaptic) knobs.
"""
import json
import sys
from pathlib import Path
import os
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'w2l_cpg'))
sys.path.insert(0,str(HERE.parent/'reports_20261016/w2l_wean'))
import _w2l_wean_eval as EV
import build_w2l_net as B

class Captured(Exception):
    pass

if __name__=='__main__':
    source=json.loads((HERE/'trial458_d30.json').read_text())
    knobs=source['job']['knobs']
    EV.T.KV['net']='bilateralrg2'
    original=B.build
    captured={}
    def spy(*args,**kwargs):
        n=original(*args,**kwargs)
        captured['g']=n.compiled.g_max_non.copy()
        captured['e']=n.compiled.del_e.copy()
        captured['i']=n.compiled.input_connectivity.copy()
        raise Captured
    B.build=spy
    def compile_via_pipeline(values):
        try:
            EV.evaluate(values,csv_path=str(HERE/'unused.csv'))
        except Captured:
            return {key:value.copy() for key,value in captured.items()}
        raise AssertionError('Did not traverse the real build path')
    baseline=compile_via_pipeline(knobs)
    results={}
    keys=list(source['expected_gains'])+['ch_ank','ch_contra','comm']
    for key in dict.fromkeys(keys):
        changed=dict(knobs)
        changed[key]=changed[key]*1.01 if changed[key] else .01
        actual=compile_via_pipeline(changed)
        differences={k:int(np.count_nonzero(actual[k]!=baseline[k])) for k in baseline}
        results[key]=differences
    # Unknown knobs must fail explicitly instead of being discarded.
    try:
        EV.evaluate(dict(knobs,nonexistent_knob=1.))
    except ValueError:
        pass
    else:
        raise AssertionError('Unknown knob silently accepted')
    inert=[key for key,d in results.items() if sum(d.values())==0]
    output=dict(conductance_dimensions=len(results),matrix_changes=results,
                inert_dimensions=inert,
                runtime_sensor_parameters=['ctn_amp','cref','ibnA'],
                unknown_knob_rejected=True)
    (HERE/'plumbing_audit.json').write_text(json.dumps(output,indent=2))
    print(len(results)-len(inert),'conductance dimensions change compiled matrices; inert:',inert)
