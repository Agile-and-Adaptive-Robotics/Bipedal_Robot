"""Bounded v1/v2 validator parity and throughput regression; no body steps."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import argparse,ast,hashlib,json,sys,timeit,zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import integrated_hunt_design as I


def run(root):
    root=Path(root);manifest=json.loads((root/'v1_frozen/MANIFEST.json').read_text());archive=Path(manifest['archive'])
    with zipfile.ZipFile(archive) as z:
        member=next(q for q in z.namelist() if q.endswith('/fsa_design_20261101/ibexc_coincidence.py'));old_bytes=z.read(member)
    expected=manifest['profile_bound_sources'][str(Path(I.IB.__file__))]
    if hashlib.sha256(old_bytes).hexdigest()!=expected:raise AssertionError('Frozen v1 source hash mismatch')
    old_source=old_bytes.decode('utf-8').replace('\r\n','\n');new_source=Path(I.IB.__file__).read_text()
    old_expression='not np.array_equal(wrapper.erange,(c.e_hi-c.e_lo)[wrapper.post,wrapper.pre])';new_expression='not np.array_equal(wrapper.erange,c.e_hi[wrapper.post,wrapper.pre]-c.e_lo[wrapper.post,wrapper.pre])'
    if old_source.replace(old_expression,new_expression)!=new_source:raise AssertionError('Change extends beyond declared validator expression')
    tree=ast.parse(old_source);function=next(q for q in tree.body if isinstance(q,ast.FunctionDef) and q.name=='verify_coincidence');env=I.IB.__dict__.copy();exec(compile(ast.Module(body=[function],type_ignores=[]),'frozen_v1_ibexc','exec'),env);old_verify=env['verify_coincidence'];new_verify=I.IB.verify_coincidence
    graph_parity={}
    for mode in ('function_equivalent','canonical_right_role_budget'):
        with np.load(root/mode/'COMPILED_GRAPH.npz') as old,np.load(root/'v2'/mode/'COMPILED_GRAPH.npz') as new:
            if set(old.files)!=set(new.files) or any(not np.array_equal(old[k],new[k]) for k in old.files):raise AssertionError('Numeric graph changed: '+mode)
            graph_parity[mode]={'arrays_checked':len(old.files),'all_arrays_exact':True}
    profile=root/'v2/function_equivalent/PROFILE.json';p=I.load_profile(profile);cfg=json.loads(Path(p['condition_file']).read_text());acts=[q.attrib['name'] for q in ET.parse(I.MODEL).getroot().find('actuator')];a=I.build(acts,options=I.condition_options(cfg),profile=p);b=I.build(acts,options=I.condition_options(cfg),profile=p);rng=np.random.default_rng(110202);ha=hashlib.sha256();hb=hashlib.sha256()
    with I.unbound_electrical(a),I.unbound_electrical(b):
        for step in range(128):
            u=np.zeros(len(a.inputs));u[a.ports['DRIVE']]=2.5*(0.,.5,1.,0.)[step//32]
            for name,index in a.ports.items():
                if name.startswith('HUNT_AF_'):u[index]=rng.uniform(0.,1.)
            a.step(u.copy());old_verify(a);b.step(u.copy());new_verify(b)
            if not np.array_equal(a.compiled.c.V,b.compiled.c.V) or not np.array_equal(a.compiled.c.b_gate,b.compiled.c.b_gate):raise AssertionError('Electrical trajectory changed')
            ha.update(a.compiled.c.V.tobytes());hb.update(b.compiled.c.V.tobytes())
    old_s=timeit.timeit(lambda:old_verify(b),number=100);new_s=timeit.timeit(lambda:new_verify(b),number=100)
    original=I.IB.verify_coincidence
    try:
        I.IB.verify_coincidence=old_verify;tick_old=timeit.timeit(lambda:I.verify(b,full=False),number=100)
        I.IB.verify_coincidence=new_verify;tick_new=timeit.timeit(lambda:I.verify(b,full=False),number=100)
    finally:I.IB.verify_coincidence=original
    negatives={};c=b.compiled
    for label,array,index in (('cached_range',c.erange,0),('cached_conductance',c.gmax,0),('protected_incoming_window',c.c.e_hi,(b.ibexc_coincidence_ledger['rows'][0]['index'],b.names.index(b.ibexc_coincidence_ledger['rows'][0]['candidate']['loads'][0]['name'])))):
        saved=float(array[index]);array[index]=saved+.001;rejected=[]
        for check in (old_verify,new_verify):
            try:check(b)
            except ValueError:rejected.append(True)
            else:rejected.append(False)
        array[index]=saved
        if rejected!=[True,True]:raise AssertionError('Diagnostic coverage changed: '+label)
        negatives[label]=rejected
    result={'frozen_v1_source_verified':True,'code_change_only_selected_index_subtraction':True,'compiled_graph_parity':graph_parity,'electrical_steps':128,'electrical_voltage_and_NaP_exact':True,'voltage_sha256_v1':ha.hexdigest(),'voltage_sha256_v2':hb.hexdigest(),'corruption_guards_unchanged':negatives,'timings_100_calls_s':{'IB_v1':old_s,'IB_v2':new_s,'actual_light_tick_verifier_v1':tick_old,'actual_light_tick_verifier_v2':tick_new},'speed_ratios':{'IB_validator':old_s/new_s,'actual_light_tick_verifier':tick_old/tick_new},'timing_scope':'Validator-only measurements on the same1033 neural graph; no whole-model speed claim','eliminated_temporary_bytes_per_call':len(c.V)**2*8,'replacement_range_vector_bytes':len(c.pre)*8,'body_steps':0,'physics_changes':False,'training':False}
    out=root/'v2/VALIDATOR_PARITY.json';out.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True);run(parser.parse_args().root)
