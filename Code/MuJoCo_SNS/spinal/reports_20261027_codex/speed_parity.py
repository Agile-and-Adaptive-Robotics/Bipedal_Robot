"""Full-trajectory sparse-backend check; physics remains authoritative."""
import json
import sys
import numpy as np
from recover_w2l import HERE, run_child

if __name__ == '__main__':
    dense = json.loads((HERE / 'trial458_d30.json').read_text())
    backend = sys.argv[1] if len(sys.argv) > 1 else 'sparse'
    name = 'trial458_' + backend + '_d30'
    sparse = run_child(dict(name=name, knobs=dense['job']['knobs'], backend=backend))
    a = np.load(HERE / 'trial458_d30.npz')
    b = np.load(HERE / (name + '.npz'))
    diffs = {k: float(np.max(np.abs(a[k] - b[k])))
             for k in ('q', 'com', 'pel', 'heelF', 'toeF', 'rgE', 'rigFz')}
    result = dict(max_absolute_differences=diffs,
                  bit_identical=all(v == 0. for v in diffs.values()),
                  dense_metrics=dense['metrics'], sparse_metrics=sparse['metrics'])
    (HERE / ('speed_parity_' + backend + '.json')).write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
