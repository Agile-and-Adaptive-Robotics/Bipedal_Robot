"""Supplement the source checkpoint with directly consumed neural/force data."""
import json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];SPINAL=HERE.parents[1]/'spinal';DEST=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\followup_source_snapshot')

def main():
 manifest=DEST/'MANIFEST.json';record=json.loads(manifest.read_text());paths=[SPINAL/n for n in ('joint_pf_weights.json','synergy_basis.npz','fsa_results/fsa_backsolve.npz','physics_fidelity/standing_solver.py','physics_fidelity/shared_force_law_functional_v3.json','functional_gait2392/curriculum_metrics.py','functional_gait2392/frontal_metrics.py','w2l_cpg/preallocated_exact_eval.py','w2l_cpg/exact_fast_eval.py')];paths+=[Path(__file__),HERE/'finalize_followup_record.py']
 rows=[]
 for source in paths:
  target=DEST/source.relative_to(REPO);data=source.read_bytes()
  if target.exists() and target.read_bytes()!=data:raise FileExistsError(target)
  target.parent.mkdir(parents=True,exist_ok=True)
  if not target.exists():target.write_bytes(data)
  rows.append(dict(source=str(source),snapshot=str(target),sha256=hashlib.sha256(data).hexdigest()))
 record['supplemental_direct_dependencies']=rows;manifest.write_text(json.dumps(record,indent=2),encoding='utf-8');print('Supplemented '+str(len(rows))+' direct dependencies')

if __name__=='__main__':main()
