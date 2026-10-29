"""Physical-motion curriculum metrics, independent of neural burst counts."""
import numpy as np
from scipy.signal import find_peaks

def assess(path,metrics,stage):
    z=np.load(path,allow_pickle=True);t=z['t'];q=z['q'];names=z['key_joints'].tolist()
    finite=bool(not metrics['nan'] and np.all(np.isfinite(q)))
    duration=8. if stage==3 else 16.
    complete=bool(metrics['t_end']>=duration-.01)
    if not finite or not complete:return dict(score=-1000.+metrics['t_end'],accepted=False,finite=finite,complete=complete)
    keep=(t>=5.)&(t<min(duration-.1,14.));tt=t[keep]
    if stage==3:
        tilt=float(metrics['bal_tilt_max']);sway=float(metrics['bal_sway']);load=float(metrics['bal_foot_load'])
        symmetry=float(metrics['bal_contact_sym']);fell=bool(metrics['bal_fell'])
        score=100.-2*tilt-200*sway-40*abs(symmetry-.5)+20*np.clip(load,0,1)-20*max(load-1.2,0)
        accepted=not fell and tilt<15. and sway<.06 and .3<symmetry<.7 and load>.5
        return dict(score=float(score-200*int(fell)),accepted=bool(accepted),finite=True,complete=True,
                    tilt_deg=tilt,sway_m=sway,normal_force_bw=load,symmetry=symmetry,fall=fell)
    peaks=[];ranges=[];hips=[];periods=[];counts=[];rompenalty=0.
    for side in ('r','l'):
        knee=q[keep,names.index('knee_angle_'+side)];hip=q[keep,names.index('hip_flexion_'+side)]
        indices=find_peaks(-knee,prominence=12.,distance=225)[0]
        peaks.append(tt[indices]);counts.append(len(indices));ranges.append(float(np.ptp(knee)));hips.append(float(np.ptp(hip)))
        periods.extend(np.diff(tt[indices]).tolist())
        rompenalty+=max(0.,ranges[-1]-100.)+max(0.,hips[-1]-80.)
    phases=[];coverage=[]
    for a,b in ((0,1),(1,0)):
        valid=0
        for begin,end in zip(peaks[a][:-1],peaks[a][1:]):
            candidates=peaks[b][(peaks[b]>begin)&(peaks[b]<end)]
            if len(candidates)==1:phases.append(float((candidates[0]-begin)/(end-begin)));valid+=1
        coverage.append(valid/max(1,len(peaks[a])-1))
    phase_error=float(np.mean(np.abs(np.array(phases)-.5))) if phases else .5
    cov=float(min(coverage));period=float(np.median(periods)) if periods else 0.
    cv=float(np.std(periods)/max(np.mean(periods),1e-6)) if periods else 2.
    motion=float(np.mean(np.minimum(ranges,60.))/60.+np.mean(np.minimum(hips,35.))/35.)
    # A continuous objective avoids a flat study of identical failure sentinels.
    score=30*min(min(counts)/3.,1.)+25*cov+35*(1.-2*phase_error)+10*motion-0.7*rompenalty
    score-=15*max(.8-period,0)+10*max(period-2.5,0)+15*min(cv,2.)
    score-=10*abs(ranges[0]-ranges[1])/max(max(ranges),1.)
    accepted=min(counts)>=3 and min(ranges)>=20 and max(ranges)<=100 and min(hips)>=10 and max(hips)<=80
    accepted=accepted and cov>=.8 and phase_error<.15 and .8<=period<=2.5 and cv<.35
    extra={}
    if stage==5:
        k=metrics.get('kine') or {};cycles=[k.get('n_cycles_r',0),k.get('n_cycles_l',0)]
        duty=[k.get('duty_r',1.),k.get('duty_l',1.)]
        # Contact alone or harness-driven translation cannot qualify.
        contact_ok=min(cycles)>=3 and all(.35<d<.85 for d in duty) and metrics['tilt_max']<30
        accepted=accepted and contact_ok
        score+=.08*max(-400.,float(metrics['kine_score']))
        extra=dict(contact_cycles=cycles,contact_duty=duty,kine_score=float(metrics['kine_score']),dx_m=float(metrics['dx']),contact_gate=contact_ok)
    return dict(score=float(score),accepted=bool(accepted),finite=True,complete=True,knee_cycles=counts,
                knee_ranges_deg=ranges,hip_ranges_deg=hips,cycle_phase_error=phase_error,
                cycle_coverage=cov,period_s=period,period_cv=cv,left_right_phases=phases,**extra)
