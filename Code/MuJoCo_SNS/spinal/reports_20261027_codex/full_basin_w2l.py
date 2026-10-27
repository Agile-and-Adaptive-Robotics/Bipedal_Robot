"""Perturb every numeric saved model/protocol setting, at a fixed 30 s horizon."""
import concurrent.futures
import json
import numpy as np
from recover_w2l import HERE,run_child

if __name__=='__main__':
    config=json.loads((HERE/'w2l_verified_config.json').read_text())
    base=config['knobs']
    rng=np.random.default_rng(20261028)
    jobs=[]
    for i in range(12):
        knobs={k:(v*rng.uniform(.99,1.01) if k!='dur_walk' else v)
               for k,v in base.items()}
        jobs.append(dict(name=f'full_basin_{i:02}',knobs=knobs,
                         backend='exact_fast',kind='all_saved_scalars_basin'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        results=list(pool.map(run_child,jobs))
    (HERE/'full_basin_results.json').write_text(json.dumps(results,indent=2))
    print('FULL BASIN',sum(r['metrics']['survive'] for r in results),'/',len(results))
