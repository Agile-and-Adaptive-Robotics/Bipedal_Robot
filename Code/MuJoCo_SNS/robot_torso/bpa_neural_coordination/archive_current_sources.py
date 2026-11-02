"""Archive lightweight code/configuration provenance, not duplicate model meshes."""
import hashlib,json,shutil,sys,socket,importlib.metadata
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];SPINAL=HERE.parents[1]/'spinal'
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 dest=ROOT/'source_snapshot';dest.mkdir(exist_ok=True);files=list(HERE.glob('*.py'))
 files += [SPINAL/n for n in ('runner.py','params.py','build_network.py','build_network_syn6.py','build_network_w2lvar.py','muscle_map.py','synergy_model.py','synergy_basis.npz','kine_ref.py')]
 files += [SPINAL/'w2l_cpg'/n for n in ('exact_fast_eval.py','preallocated_exact_eval.py','sparse_pairwise_eval.py')]
 files += [SPINAL/'physics_fidelity'/n for n in ('standing_solver.py','shared_force_law_functional_v3.json')]
 files += [SPINAL/'functional_gait2392'/n for n in ('curriculum_metrics.py','frontal_metrics.py')]
 files += [HERE.parents[1]/'bpa_muscle.py',HERE/'GAIT_SYN6_STANDING_MANIFEST.json']
 files += [ROOT/'rustuna_syn6_v1/BEST.json',ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json',ROOT/'virtual_bpa_design/HELDOUT_SLOW_WEAK_MEMORY.json',ROOT/'virtual_bpa_design/HELDOUT_MAPERROR.json',ROOT/'virtual_bpa_design/HELDOUT_DELAY.json',ROOT/'pinned_preparation/PREPARATION_REFERENCE.json',ROOT/'pinned_preparation/SELECTED_RUNG.json']
 entries=[];missing=[]
 for source in files:
  if not source.exists():missing.append(str(source));continue
  relative=source.relative_to(REPO) if source.is_relative_to(REPO) else Path('configuration')/source.relative_to(ROOT)
  target=dest/relative;target.parent.mkdir(parents=True,exist_ok=True)
  if target.exists() and sha(target)!=sha(source):
   if source.resolve()==Path(__file__).resolve():shutil.copyfile(source,target)
   else:raise FileExistsError('Preserve archived source version: '+str(target))
  if not target.exists():shutil.copyfile(source,target)
  entries.append(dict(source=str(source),snapshot=str(target),sha256=sha(source)))
 versions={}
 for name in ('numpy','mujoco','sns-toolbox','numba','scipy','matplotlib'):
  try:versions[name]=importlib.metadata.version(name)
  except importlib.metadata.PackageNotFoundError:versions[name]='not found'
 native=ROOT/'pinned_preparation/syn6_prepared_nominal_v1.npz'
 with np.load(native) as z:q0=z['qfull'][0].copy()
 np.savez_compressed(dest/'INITIAL_QPOS.npz',qpos=q0)
 model=REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml'
 manifest=dict(machine=socket.gethostname(),python_executable=sys.executable,python_version=sys.version,packages=versions,files=entries,missing_optional_files=missing,generated_model_dependency=dict(path=str(model),sha256=sha(model)),scope='Lightweight provenance archive only. Model XML/meshes and external package binaries are not copied; matching runtime dependencies are still required. Initial qpos is preserved separately. Existing snapshots refuse changed-source overwrite.')
 (dest/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps(dict(snapshot=str(dest),files=len(entries),missing=missing,packages=versions),indent=2))

if __name__=='__main__':main()
