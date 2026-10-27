import concurrent.futures
import json
from recover_w2l import HERE,run_child

if __name__=='__main__':
    original=json.loads((HERE/'trial458_d30.json').read_text())
    audit=json.loads((HERE/'plumbing_audit.json').read_text())
    knobs=dict(original['job']['knobs'])
    inert={k:knobs.pop(k) for k in audit['inert_dimensions']}
    (HERE/'w2l_verified_config.json').write_text(json.dumps(dict(
        knobs=knobs,backend='exact_fast',
        effective_conductance_dimensions=56,runtime_sensor_dimensions=3,
        excluded_legacy_inert_keys=inert,
        scope='Harness-supported W2L stepping; contact-event gate, not physiological walking'),indent=2))
    jobs=[dict(name='canonical_repeat',knobs=knobs,backend='exact_fast')]
    # Test the more balanced contact-count candidate as well as the score winner.
    candidate=json.loads((HERE/'trial360_d30.json').read_text())['job']['knobs']
    import numpy as np
    rng=np.random.default_rng(20261027)
    live=[k for k in audit['matrix_changes'] if k not in audit['inert_dimensions']]
    live+=audit['runtime_sensor_parameters']
    for i in range(12):
        changed=dict(candidate)
        for k in live:
            changed[k]*=rng.uniform(.99,1.01)
        jobs.append(dict(name=f'balanced_basin_{i:02}',knobs=changed,
                         backend='exact_fast',kind='balanced_candidate_basin'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results=list(pool.map(run_child,jobs))
    (HERE/'balanced_basin_results.json').write_text(json.dumps(results,indent=2))
    a=np.load(HERE/'trial458_d30.npz'); b=np.load(HERE/'canonical_repeat.npz')
    assert all(np.array_equal(a[k],b[k]) for k in
               ('q','com','pel','heelF','toeF','rgE','rigFz'))
    print('Canonical config exactly reproduces the reference trajectory.')
