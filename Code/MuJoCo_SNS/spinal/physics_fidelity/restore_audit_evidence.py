"""Reproduce lost volatile scratch evidence in the durable ignored audit folder.

Does not export or overwrite any curated force-law parameter JSON.
"""
import sys,json,subprocess,hashlib,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
HERE=Path(__file__).resolve().parent;OUT=HERE/'runs/audit_20261030';OUT.mkdir(parents=True,exist_ok=True)
names=('benchmark.py','contact_benchmark.py','opensim_static_check.py','follower_inertia_audit.py','bilateral_moment_arm_audit.py','conversion_cache_static_check.py','conversion_activation_check.py','shared_force_surrogate_audit.py')
def run(name):
 with (OUT/(name+'.log')).open('w') as log:code=subprocess.run([sys.executable,str(HERE/name)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,code
before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('shared_force_law_functional_v*.json')}
with ThreadPoolExecutor(max_workers=2) as pool:rows=dict(pool.map(run,names))
after={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('shared_force_law_functional_v*.json')}
assert before==after
manifest=dict(reason='Earlier tmp audit files disappeared; cause unknown. These are regenerated evidence, not recovery of identical historical raw-file hashes.',exit_codes=rows,curated_parameter_hashes=after)
(OUT/'REGENERATED_EVIDENCE.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest),flush=True)
assert not any(rows.values()),rows
