"""Run one explicit joint-integration probe with the existing fast evaluator."""
import os, sys, json, runpy, hashlib
from pathlib import Path
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parent), str(H.parent / 'w2l_cpg')]
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
os.environ['AARL_NET'] = 'w2lvar'
config = Path(sys.argv[1]).resolve()
sys.argv[1] = str(config)
cfg = json.loads(config.read_text())
if cfg.get('network_variant') != 'deng2023_joint_integration':
    raise ValueError('Configuration must explicitly select the experimental variant.')
import network as F
import deng2023_joint_pf as D
capture = {}
def build(*args, **kwargs):
    kwargs['options'] = {**kwargs.get('options', {}), 'joint_mix_gain': cfg.get('joint_mix_gain', 1.)}
    net = D.build(*args, **kwargs)
    capture.update(joint_integration_edges=net.joint_integration_edges,
                   bi_integrators=net.bi_integrators)
    return net
F.build = build
runpy.run_path(str(H / 'curriculum_eval.py'), run_name='__main__')
resultpath = config.with_suffix('.result.json')
result = json.loads(resultpath.read_text())
result.update(network_variant=cfg['network_variant'], joint_mix_gain=cfg.get('joint_mix_gain', 1.), **capture)
result['source_hashes'].update({p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in (Path(__file__), H / 'deng2023_joint_pf.py')})
result['interpretation'] = 'Untuned human joint-integration probe; no Deng replication or curriculum completion claim.'
resultpath.write_text(json.dumps(result, indent=2))
