"""Stage-appropriate physical and neural gates; no sentinel-score winners."""
import os,sys,json
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from scipy.signal import find_peaks
OUT=Path(__file__).resolve().parent.parent/'reports_20261028_functional_pf'
results=[]
for path in sorted((OUT/'trials').glob('*.json')):
    if path.stem.endswith('_assessment'):continue
    summary=json.loads(path.read_text());npz=path.with_suffix('.npz')
    if not npz.exists():continue
    d=np.load(npz,allow_pickle=True);keep=(d['t']>=2)&(d['t']<13);t=d['t'][keep]
    names=d['neuro_names'].tolist();joints=d['key_joints'].tolist()
    bursts=[];excursions=[];cycles=[];diff=[]
    for side in ('r','l'):
        signal=d['neuro'][keep,names.index('RG_F_'+side)]-d['neuro'][keep,names.index('RG_E_'+side)]
        diff.append(signal);bursts.append(int(len(find_peaks(signal,prominence=1.,distance=150)[0])))
        knee=d['q'][keep,joints.index('knee_angle_'+side)]
        excursions.append(float(np.ptp(knee)));cycles.append(int(len(find_peaks(-knee,prominence=5.,distance=150)[0])))
    corr=float(np.corrcoef(diff)[0,1]) if min(np.std(x) for x in diff)>1e-9 else 0.
    gates=dict(finite=not summary['metrics']['nan'],completed=summary['metrics']['t_end']>=(7.9 if summary['stage']=='stand' else 15.9),
        bilateral_rg_rhythm=min(bursts)>=3,antiphase=corr<-.5,bilateral_knee_motion=min(excursions)>=15.,bilateral_physical_cycles=min(cycles)>=3)
    late=(d['t']>=5)&(d['t']<13);lt=d['t'][late];apices=[];late_ranges=[]
    for side in ('r','l'):
        knee=d['q'][late,joints.index('knee_angle_'+side)]
        late_ranges.append(float(np.ptp(knee)));apices.append(lt[find_peaks(-knee,prominence=5.,distance=175)[0]])
    phases=[]
    for begin,end in zip(apices[0][:-1],apices[0][1:]):
        candidates=apices[1][(apices[1]>begin)&(apices[1]<end)]
        if len(candidates)==1:phases.append(float((candidates[0]-begin)/(end-begin)))
    coverage=len(phases)/max(1,len(apices[0])-1)
    gates['bounded_late_knee_excursion']=max(late_ranges)<=110.
    gates['alternating_physical_knee_cycles']=len(phases)>=2 and coverage>=.8 and all(.2<p<.8 for p in phases)
    result=dict(run=path.stem,stage=summary['stage'],variant=summary['variant'],gates=gates,rg_bursts=bursts,
        knee_excursions_deg=excursions,knee_cycle_counts=cycles,rg_correlation=corr,
        late_knee_excursion_deg=late_ranges,left_knee_apex_phases=phases,matched_knee_cycle_fraction=coverage,
        kine_score=summary['metrics']['kine_score'],air_pattern_pass=all(gates.values()) if summary['stage'].startswith('air') else None,
        scope='Air-pattern diagnostics only; ground and standing outcomes require their own load, posture, contact, and progression gates')
    results.append(result)
(OUT/'stage_assessments.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
