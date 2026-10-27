"""Frozen full-conductance stimulus and deletion experiment."""
import concurrent.futures
import json
from recover_w2l import HERE, run_child


def main():
    item = json.loads((HERE / 'trial458_d30_ws25.json').read_text())
    base = item['job']['knobs']
    conditions = dict(control={},
        stim_mild=dict(stimN=1.5, stimT=8., stimDur=1.),
        stim_strong=dict(stimN=3., stimT=8., stimDur=1.),
        push_fwd=dict(pushN=40., pushT=8., pushAxis='+x'),
        push_back=dict(pushN=40., pushT=8., pushAxis='-x'),
        del_crossed_heel=dict(ch_contra=0.),
        del_Ia={k: 0. for k in base if k.startswith('ia_')},
        del_Renshaw={k: 0. for k in base if k.startswith('rc_')},
        del_heel={k: 0. for k in base if k.startswith('heel_')},
        del_commissural=dict(comm=0.),
        del_Ib=dict({k: 0. for k in base if k.startswith('ib_')}, ibnA=0.))
    conditions['del_all_feedback'] = dict(
        {k: 0. for k in base if k.startswith(('ia_', 'ii_', 'ib_', 'heel_'))},
        ctn_amp=0., ibnA=0., ch_ank=0., ch_contra=0.)
    jobs = [dict(name='panel_' + name, knobs=dict(base, **over),
                 condition=name, kind='frozen_perturbation')
            for name, over in conditions.items()]
    # Preserve the protocol's text-valued push-axis setting separately;
    # the historical evaluator accepts numeric knobs only.
    for job in jobs:
        if 'pushAxis' in job['knobs']:
            job['text_knobs'] = dict(pushAxis=job['knobs'].pop('pushAxis'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(run_child, jobs))
    (HERE / 'perturbation_results.json').write_text(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
