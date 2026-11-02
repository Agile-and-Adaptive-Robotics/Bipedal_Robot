"""Read-only numeric comparisons and measured timing ratios for saved runs."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np


def compare_npz(a,b):
    with np.load(a,allow_pickle=False)as x,np.load(b,allow_pickle=False)as y:
        rows=[];skipped=[]
        for name in sorted(set(x.files)|set(y.files)):
            if name not in x.files or name not in y.files:
                rows.append(dict(field=name,present_both=False,exact=False));continue
            try:u=x[name];v=y[name]
            except ValueError:
                skipped.append(dict(field=name,reason='Object/pickled metadata, not a numeric array'));continue
            if u.dtype.kind not in'biufc'or v.dtype.kind not in'biufc':
                skipped.append(dict(field=name,reason='Text/label metadata'));continue
            row=dict(field=name,present_both=True,shape_a=list(u.shape),shape_b=list(v.shape),
                     dtype_a=str(u.dtype),dtype_b=str(v.dtype),same_shape=u.shape==v.shape,same_dtype=u.dtype==v.dtype)
            if u.shape!=v.shape:row.update(exact=False,max_abs_error=None)
            else:
                finite=np.isfinite(u)&np.isfinite(v)
                same_masks=bool(np.array_equal(np.isnan(u),np.isnan(v))and np.array_equal(np.isposinf(u),np.isposinf(v))and np.array_equal(np.isneginf(u),np.isneginf(v)))
                row.update(exact=bool(np.array_equal(u,v,equal_nan=True)),nonfinite_masks_exact=same_masks,
                           finite_entries=int(np.count_nonzero(finite)),
                           max_abs_error=float(np.max(abs(u[finite].astype(float)-v[finite].astype(float))))if np.any(finite)else None)
            rows.append(row)
    return dict(file_a=str(a),file_b=str(b),file_sha256_a=hashlib.sha256(a.read_bytes()).hexdigest(),
                file_sha256_b=hashlib.sha256(b.read_bytes()).hexdigest(),
                numeric_fields=rows,skipped_metadata=skipped,
                all_numeric_fields_exact=bool(rows)and all(r.get('exact',False)for r in rows),
                interpretation='Exact numeric values with equal NaN masks; all-NaN unlogged stimulus arrays remain unlogged, not evidence of zero applied stimulus.')


def timing(root,label):
    record=root/'timings'/label/'RESULT.json';result=root/'configurations'/(label+'.result.json')
    a=json.loads(record.read_text());b=json.loads(result.read_text())
    return dict(label=label,timing_file=str(record),result_file=str(result),
                whole_process_wall_s=a['total_process_wall_s'],inner_evaluator_wall_s=b['wall_s'],
                accumulated_neural_call_wall_s=b['compiled'].get('neural_step_wall_s'),
                reported_protocol_s=b.get('simulated_s'),backend=b.get('backend'),
                accepted=b['assessment']['accepted'],exit_code=a['exit_code'],timed_out=a['timed_out'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root.resolve()
    definitions=[('backend','source_baseline_s1','source_native_s1'),('validator','fsa_equivalent_s2','fsa_v2_equivalent_s2')];pairs=[]
    suffixes=('.npz','.assistance.npz','.balance.npz','.contact_events.npz','.regional_contacts.npz','.mechanics.npz','.neurons.npz')
    for kind,left,right in definitions:
        files=[compare_npz(root/'figures'/left/'aligned_trace.npz',root/'figures'/right/'aligned_trace.npz')]
        for suffix in suffixes:
            x=root/'configurations'/(left+suffix);y=root/'configurations'/(right+suffix)
            if x.exists()and y.exists():files.append(compare_npz(x,y))
            elif x.exists()or y.exists():files.append(dict(file_a=str(x),file_b=str(y),missing_counterpart=True,all_numeric_fields_exact=False))
        lt=timing(root,left);rt=timing(root,right)
        ratios={key:rt[key]/lt[key]for key in('whole_process_wall_s','inner_evaluator_wall_s','accumulated_neural_call_wall_s')if lt.get(key)and rt.get(key)}
        if kind=='validator':ratios={key:1/value for key,value in ratios.items()}
        pairs.append(dict(kind=kind,run_a=left,run_b=right,files=files,
                          all_saved_numeric_arrays_exact=all(f['all_numeric_fields_exact']for f in files),
                          timings=[lt,rt],observed_speedup=ratios,
                          ratio_direction='native/accelerated'if kind=='backend'else'validator_v1/validator_v2'))
    result=dict(machine='easteregg2',pairs=pairs,comparison_scope='Same saved verification cohorts; aligned capture and every available matching numeric NPZ sidecar. Text, labels, object metadata and result provenance are excluded from numeric equality.',
                timing_scope='Observed paired fresh-interpreter runs with capture and diagnostics; whole process includes startup/compilation/evaluation/assessment/output. Accumulated neural calls are not a pure-kernel benchmark. No universal speedup, statistical distribution, or faster-than-real-time claim.',
                no_new_simulation=True,all_pairs_exact=all(pair['all_saved_numeric_arrays_exact']for pair in pairs))
    output=root/'figures';output.mkdir(exist_ok=True)
    (output/'BACKEND_AND_VALIDATOR_COMPARISON.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    lines=['# Saved backend and validator comparison','','Machine: **easteregg2**. Existing saved runs only; no new simulation.','',
           '## Numeric parity','', '| Pair | Numeric array fields | Files | Exact numeric parity |','| --- | ---: | ---: | --- |']
    for pair in pairs:
        count=sum(len(f.get('numeric_fields',[]))for f in pair['files']);lines.append(f"| {pair['run_a']} / {pair['run_b']} | {count} | {len(pair['files'])} | {pair['all_saved_numeric_arrays_exact']} |")
    lines+=['','Equality includes times, selected actual neural voltages, qpos/qvel, act/ctrl, and all recorded force arrays, plus numeric source/diagnostic sidecars. NaN masks must agree. Unlogged all-NaN stimuli are not treated as zero stimuli. Labels and metadata are not used to establish numerical equality. The JSON gives every field’s shape, dtype, finite-entry count, maximum absolute difference, and equality result.','',
            '## Measured wall time','', '| Run | Whole process (s) | Inner evaluator (s) | Neural calls (s) |','| --- | ---: | ---: | ---: |']
    for pair in pairs:
        for t in pair['timings']:lines.append(f"| {t['label']} | {t['whole_process_wall_s']:.6f} | {t['inner_evaluator_wall_s']:.6f} | {t['accumulated_neural_call_wall_s']:.6f} |")
    for pair in pairs:
        r=pair['observed_speedup'];lines+=['',f"{pair['ratio_direction']}: whole process **{r['whole_process_wall_s']:.3f}×**, inner evaluator **{r['inner_evaluator_wall_s']:.3f}×**, accumulated neural calls **{r['accumulated_neural_call_wall_s']:.3f}×**."]
    lines+=['','These are measured paired verification runs, including their capture/diagnostic costs and fresh-process startup. They are not a universal 200× claim or a pure neural-kernel benchmark. The two validation versions retain their own configuration/source lineage; equality concerns actual saved numeric traces. Protocol duration and native physics clock are kept distinct because held-pose preparation does not advance physical time. No curriculum acceptance is changed by this comparison.','',
            '## Per-field record','']
    for pair in pairs:
        lines+=[f"### {pair['run_a']} / {pair['run_b']}",'']
        for f in pair['files']:
            lines+=[f"File: `{f['file_a']}` versus `{f['file_b']}`",'', '| Field | Shape A / B | Maximum absolute error | Exact |','| --- | --- | ---: | --- |']
            for r in f.get('numeric_fields',[]):lines.append(f"| {r['field']} | {r.get('shape_a')} / {r.get('shape_b')} | {r.get('max_abs_error')} | {r['exact']} |")
            lines+=['']
    (output/'BACKEND_AND_VALIDATOR_COMPARISON.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(dict(all_pairs_exact=result['all_pairs_exact'],pairs=[dict(kind=q['kind'],exact=q['all_saved_numeric_arrays_exact'],speedup=q['observed_speedup'])for q in pairs]),indent=2))


if __name__=='__main__':main()
