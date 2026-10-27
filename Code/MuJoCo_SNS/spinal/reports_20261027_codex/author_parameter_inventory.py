"""Inventory the supplied optimized parameters without collapsing edges."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
path=REPO/'tmp/codex_dirusso_original/import/initializations/spinal.par'
parameters={line.split()[0]:float(line.split()[1]) for line in path.read_text().splitlines() if line.strip()}
categories={
    'frequency':[k for k in parameters if k=='phase_oscillator.omega'],
    'initial_state_offsets':[k for k in parameters if k.endswith('.offset') and not k.startswith('phase_')],
    'initial_phase_offsets':[k for k in parameters if k in ('phase_left.offset','phase_right.offset')],
    'balance':[k for k in parameters if k.startswith('S00011.')],
    'pattern_weights':[k for k in parameters if k.startswith(tuple(f'P{i}.' for i in range(5)))],
    'reflex_connections':[k for k in parameters if '->' in k and not k.startswith(tuple(f'P{i}.' for i in range(5)))],
    'sensory_length_offsets':[k for k in parameters if k.endswith('SN_II.w0')],
}
claimed=[k for keys in categories.values() for k in keys]
assert len(claimed)==len(set(claimed))==len(parameters)
out=dict(total=len(parameters),counts={k:len(v) for k,v in categories.items()},
          categories=categories,parameters=parameters)
(HERE/'author_parameter_inventory.json').write_text(json.dumps(out,indent=2))
print(out['total'],out['counts'])
