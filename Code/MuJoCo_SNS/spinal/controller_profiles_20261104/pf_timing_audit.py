"""Audit actual PF transmission timing before allowing MN readout fitting."""
import argparse,json,hashlib,sys
from pathlib import Path
import numpy as np
from functools import lru_cache
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))

@lru_cache(maxsize=1)
def left_phase_offset():
 import kine_ref
 return float(kine_ref.load_reference()['lag_rl'])

def primitive_target(fsa,cfg,side,local_time,k):
 # The current composed syn6 profile mirrors RIGHT PF circuit/readout onto
 # the left. Independently factorized original left NMF columns are not that
 # mirrored target contract and must not be treated as matching S labels.
 shift=left_phase_offset()*cfg['target_period_s']if side=='l'else 0.
 query=cfg['target_start_s']+(np.asarray(local_time)-shift)%cfg['target_period_s']
 return np.clip(np.interp(query,fsa['time'],fsa['r_source'][:,k]),0.,1.)

def features(y):
 active=y>.2;on=np.flatnonzero(active&~np.roll(active,1));off=np.flatnonzero(~active&np.roll(active,1))
 phase=np.arange(len(y))/len(y);coefficient=np.mean((y-y.mean())*np.exp(-2j*np.pi*phase))
 return dict(minimum=float(y.min()),maximum=float(y.max()),excursion=float(np.ptp(y)),duty_fraction=float(active.mean()),onset_phase_percent=(on/len(y)*100).tolist(),offset_phase_percent=(off/len(y)*100).tolist(),peak_phase_percent=float(100*phase[np.argmax(y)]),phase_defined=bool(np.ptp(y)>.05),fundamental_phase_deg=float(np.angle(coefficient,deg=True)))

def effective_rank(values):
 centered=values-values.mean(axis=0);s=np.linalg.svd(centered,compute_uv=False);energy=s*s
 return dict(singular_values=s.tolist(),rank_for99percent_variance=int(np.searchsorted(np.cumsum(energy)/max(energy.sum(),1e-30),.99)+1))

def main():
 p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 dm=json.loads((a.fit/'MANIFEST.json').read_text());cm=json.loads((Path(dm['fit'])/'MANIFEST.json').read_text());cfg=cm['configuration'];family=cfg['family'];T=cfg['target_period_s'];z=np.load(a.fit/'BEST_TRACE.npz');names=z['pf_names'].tolist();phase=np.arange(200)/200.;sample_t=phase*T;rows=[];curves={};collapse=[];rank={}
 fsa_path=HERE.parent/'fsa_results/fsa_backsolve.npz';fsa=np.load(fsa_path)
 for side in('r','l'):
  selected=[(i,n)for i,n in enumerate(names)if n.endswith('_'+side)];actual=[];targets=[]
  for col,name in selected:
   cycles=np.array([np.interp(sample_t+n*T,z['t'],np.clip(z['pf'][:,col]/5.,0.,1.))for n in range(2,11)])
   y=cycles.mean(axis=0);row=dict(cell=name,actual=features(y),cycle_duty_mean=float(np.mean((cycles>.2).mean(axis=1))),cycle_duty_SD=float(np.std((cycles>.2).mean(axis=1),ddof=1)))
   target=None
   if family=='syn6':
    k=int(name.split('_')[1][1:])-1;target=primitive_target(fsa,cfg,side,sample_t,k);wanted=features(target)
    lag=(wanted['fundamental_phase_deg']-row['actual']['fundamental_phase_deg']+180.)%360.-180.
    row.update(target=wanted,duty_error_fraction=row['actual']['duty_fraction']-wanted['duty_fraction'],fundamental_lag_deg=float(lag),peak_phase_error_deg=float((3.6*(row['actual']['peak_phase_percent']-wanted['peak_phase_percent'])+180.)%360.-180.),mean_absolute_activation_error=float(np.mean(abs(y-target))),maximum_absolute_activation_error=float(np.max(abs(y-target))),onset_offset_count_match=len(wanted['onset_phase_percent'])==len(row['actual']['onset_phase_percent'])and len(wanted['offset_phase_percent'])==len(row['actual']['offset_phase_percent']))
    # Provisional, declared calibration gate. Never silently promote an MN fit.
    row['timing_gate']=bool(abs(lag)<=30. and abs(row['duty_error_fraction'])<=.15 and row['onset_offset_count_match']and row['actual']['excursion']>=.2)
    targets.append(target)
   curves[name]=(y,target);actual.append(y);rows.append(row)
  matrix=np.array(actual).T;rank[side]=dict(actual=effective_rank(matrix));corr=np.corrcoef(matrix.T)
  if targets:rank[side]['target']=effective_rank(np.array(targets).T);target_corr=np.corrcoef(np.array(targets))
  for i in range(len(selected)):
   for j in range(i+1,len(selected)):
    if corr[i,j]>.98:collapse.append(dict(cells=[selected[i][1],selected[j][1]],actual_correlation=float(corr[i,j]),target_correlation=float(target_corr[i,j])if targets else None))
 report=dict(family=family,fit=str(a.fit),reference_period_s=T,clock='Original calibration time; no channel-specific phase shift and no heel/toe display time warp',target_contract='Composed syn6 mirrors right basis. Left uses right primitive with measured contralateral loading phase offset, not the independently factorized original left NMF columns.',activation='clip(PF voltage /5mV,0,1); active threshold0.2',rows=rows,highly_correlated_pairs=collapse,effective_rank=rank,PF_timing_verified=bool(all(r.get('timing_gate',False)for r in rows)),gate_scope='syn6 provisional engineering gate: |fundamental lag|<=30deg, |duty error|<=0.15, matching burst counts and excursion>=0.2. s3k has no declared numerical channel profiles; cannot certify timing from labels alone.',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in(Path(__file__),fsa_path,a.fit/'BEST_TRACE.npz')})
 (a.output/'PF_TIMING_AUDIT.json').write_text(json.dumps(report,indent=2))
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'Arial','font.size':10})
 n=len(rows)//2;fig,axes=plt.subplots(n,2,figsize=(11,2*n),sharex=True,sharey=True,constrained_layout=True)
 for col,side in enumerate(('r','l')):
  subset=[r for r in rows if r['cell'].endswith('_'+side)]
  for ax,row in zip(axes[:,col],subset):
   actual,target=curves[row['cell']];ax.plot(phase*100,actual,color='#4477AA',label='Actual PF transmission')
   if target is not None:ax.plot(phase*100,target,'--',color='#222222',label='Intended primitive')
   ax.axhline(.2,color='#BBBBBB',lw=.6);ax.set_title(row['cell']);ax.set_ylim(-.05,1.1);ax.grid(alpha=.15);ax.set_ylabel('Activation fraction')
  axes[-1,col].set_xlabel('Original calibration cycle (%)')
 axes[0,0].legend(fontsize=8);fig.suptitle(f'{family}: PF timing audit before MN readout fitting\nSame reference clock; no individual phase realignment')
 for ext in('png','pdf','svg'):fig.savefig(a.output/('PF_TIMING_AUDIT.'+ext),dpi=170)
 plt.close(fig)
 print(json.dumps(dict(family=family,verified=report['PF_timing_verified'],rank=rank,channels=[dict(cell=r['cell'],duty=r['actual']['duty_fraction'],target_duty=r.get('target',{}).get('duty_fraction'),lag=r.get('fundamental_lag_deg'),passed=r.get('timing_gate'))for r in rows])))
if __name__=='__main__':main()
