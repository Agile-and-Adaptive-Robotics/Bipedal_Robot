"""Topology, opt-in isolation, and exact-backend checks for joint integration."""
import os, sys, json, hashlib
from pathlib import Path
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parent), str(H.parent / 'w2l_cpg')]
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
os.environ['AARL_NET'] = 'w2lvar'
import numpy as np
import mujoco
import runner as R
import network as F
import deng2023_joint_pf as D
from preallocated_exact_eval import PreallocatedExactForward
from sparse_pairwise_eval import SparsePairwiseForward
model = mujoco.MjModel.from_xml_path(str(R.MODEL))
acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i) for i in range(model.nu)]
options = dict(axes=True, biarticular=True)
before = F.build(acts, options=options)
candidate = D.build(acts, options=options)
second = D.build(acts, options=options)
after = F.build(acts, options=options)
names = [p['name'] for p in candidate.net.populations]
for key, target in candidate.bi_integrators.items():
    family, side = key.rsplit('_', 1)
    incoming = [c for c in candidate.net.connections if names[c['destination']] == target]
    assert {names[c['source']] for c in incoming} == {
        template.format(side=side) for template, _ in D.JOINT_INPUTS[family]}
    assert len(incoming) == 2
    for act, info in candidate.routing['actuator_mapping'].items():
        if info['family'] == family and act.endswith('_' + side):
            pf = [names[c['source']] for c in candidate.net.connections
                  if names[c['destination']] == 'MN_' + act and names[c['source']].startswith('PF_')]
            assert target in pf
            assert not any(p.startswith(('PF_HIP-', 'PF_KNEE-', 'PF_FUNC_BI_')) for p in pf)
assert len(candidate.bi_integrators) == 6
assert len(before.compiled.V) == len(after.compiled.V) == 980
assert len(candidate.compiled.V) == 962
dense = PreallocatedExactForward(candidate.compiled)
sparse = SparsePairwiseForward(second.compiled)
rng = np.random.default_rng(20231028)
for step in range(1000):
    u = rng.uniform(0., 1.5, len(candidate.inputs))
    dense.forward(u); sparse.forward(u)
    assert np.array_equal(dense.V, sparse.V), step
    assert np.array_equal(dense.c.b_gate, sparse.c.b_gate), step
    before.compiled.forward(u); after.compiled.forward(u)
    assert np.array_equal(before.compiled.V, after.compiled.V), step
assert np.all(np.isfinite(sparse.V))
out = H.parent / 'reports_20261028_functional_pf/literature_training_review'
result = dict(neurons=962, inputs=len(candidate.inputs), joint_integrators=candidate.bi_integrators,
    source_edges=candidate.joint_integration_edges, verification_steps=1000,
    exact_sparse_vs_preallocated=True, exact_existing_variant_before_after=True,
    native_outputs_finite=True, source_sha256=hashlib.sha256((H / 'deng2023_joint_pf.py').read_bytes()).hexdigest())
(out / 'joint_integration_verification.json').write_text(json.dumps(result, indent=2))
print(json.dumps({k: v for k, v in result.items() if k not in ('source_edges', 'joint_integrators')}, indent=2))
