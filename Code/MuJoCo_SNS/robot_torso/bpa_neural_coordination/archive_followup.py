"""Immutable lightweight follow-up source archive; mesh/model dependencies stay in place."""
import sys,json,socket,hashlib,shutil,importlib.metadata
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];SPINAL=HERE.parents[1]/'spinal';ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031');DEST=ROOT/'gait_bpa/followup_source_snapshot'

def main():
 DEST.mkdir(exist_ok=True);files=list(HERE.glob('*.py'))
 files += [SPINAL/n for n in ('runner.py','params.py','build_network.py','build_network_w2lvar.py','build_network_syn6.py','muscle_map.py','kine_ref.py')]
 files += [SPINAL/'virtual_walker_20261031'/n for n in ('spiking_sns_inputs.py','spiking_sensory.py','rustuna_transfer.py','rustuna_myo_worker.py','rustuna_runtime.py')]
 files += [SPINAL/'w2l_cpg/sparse_pairwise_eval.py',HERE.parent/'dynamic_routes.py',HERE.parent/'robotbody_routes.json']
 files += [ROOT/'gait_bpa'/n for n in ('GAIT_S3K_STANDING_MANIFEST.json','GAIT_W2LVAR_STANDING_MANIFEST.json','rustuna_s3k_v1/BEST.json','rustuna_w2lvar_v1/BEST.json','rustuna_syn6_v1/BEST.json','virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json','virtual_bpa_design/HELDOUT_SLOW_WEAK_MEMORY.json','pinned_preparation/PREPARATION_REFERENCE.json')]
 rows=[]
 for source in files:
  relative=source.relative_to(REPO) if source.is_relative_to(REPO) else Path('configuration')/source.relative_to(ROOT);target=DEST/relative;content=source.read_bytes()
  if target.exists() and target.read_bytes()!=content:raise FileExistsError('Preserve archived source: '+str(target))
  target.parent.mkdir(parents=True,exist_ok=True)
  if not target.exists():target.write_bytes(content)
  rows.append(dict(source=str(source),snapshot=str(target),sha256=hashlib.sha256(content).hexdigest()))
 versions={n:importlib.metadata.version(n) for n in ('numpy','scipy','mujoco','sns-toolbox','numba','matplotlib')}
 modelpaths=[REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody_cvt3.xml',HERE.parent/'ah001_neural.xml']
 report=dict(machine=socket.gethostname(),python_executable=sys.executable,python_version=sys.version,packages=versions,files=rows,model_dependencies=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in modelpaths],scope='Current post-follow-up source checkpoint, not a retroactive assertion that every earlier logging-only revision was identical. Numerical result files, per-run hashes, and state audits remain authoritative. Model meshes and package binaries are not duplicated; this source archive alone is not a standalone installation.')
 (DEST/'MANIFEST.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(files=len(rows),manifest=str(DEST/'MANIFEST.json')),indent=2))

if __name__=='__main__':main()
