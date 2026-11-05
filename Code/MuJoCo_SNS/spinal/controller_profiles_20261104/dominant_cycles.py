"""Audit physical principal cycles without counting secondary knee humps twice.

The historical assessor is retained. This separate check estimates period from
joint autocorrelation, then selects separated flexion peaks. It reports every
legacy secondary peak and rejects fast oscillation, frozen, or synchronous legs.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.signal import correlate,find_peaks,periodogram

def measure(t,q,hip):
    dt=float(np.median(np.diff(t)));peaks=[];periods=[];correlations=[];secondary=[];spectral=[]
    for j in range(2):
        y=q[:,j]-np.mean(q[:,j]);variance=float(np.mean(y*y))
        if variance<1e-8:return dict(accepted=False,reason='frozen_joint')
        frequency,power=periodogram(y,fs=1/dt,window='hann');dominant=float(frequency[np.argmax(power[1:])+1]);spectral.append(dominant)
        if not .35<=dominant<=1.35:return dict(accepted=False,reason='dominant_frequency_outside_walking_window',dominant_Hz=spectral)
        ac=correlate(y,y,mode='full',method='fft')[len(y)-1:]/(np.arange(len(y),0,-1)*variance)
        candidates=find_peaks(ac)[0];candidates=candidates[(candidates*dt>=.8)&(candidates*dt<=2.5)]
        if not len(candidates):return dict(accepted=False,reason='no_principal_autocorrelation_peak')
        k=int(candidates[np.argmax(ac[candidates])]);T=k*dt
        coefficient=float(np.corrcoef(y[:-k],y[k:])[0,1]);correlations.append(coefficient)
        if coefficient<.6:return dict(accepted=False,reason='nonperiodic_joint',correlation=coefficient)
        raw=find_peaks(q[:,j],prominence=12.,distance=max(1,round(.45/dt)))[0]
        picked=find_peaks(q[:,j],prominence=12.,distance=max(1,round(.65*T/dt)))[0]
        secondary.append(dict(all_prominent_peaks=len(raw),principal_peaks=len(picked)))
        peaks.append(t[picked]);periods.extend(np.diff(t[picked]).tolist())
    phases=[];coverage=[]
    for a,b in ((0,1),(1,0)):
        valid=0
        for begin,end in zip(peaks[a][:-1],peaks[a][1:]):
            p=peaks[b][(peaks[b]>begin)&(peaks[b]<end)]
            if len(p)==1:phases.append(float((p[0]-begin)/(end-begin)));valid+=1
        coverage.append(valid/max(1,len(peaks[a])-1))
    phase_error=float(np.mean(abs(np.array(phases)-.5)))if phases else .5
    T=float(np.median(periods))if periods else 0.;cv=float(np.std(periods)/max(np.mean(periods),1e-9))if periods else 99.
    knee_rom=np.ptp(q,axis=0);hip_rom=np.ptp(hip,axis=0)
    accepted=min(map(len,peaks))>=3 and min(knee_rom)>=20 and max(knee_rom)<=100 and min(hip_rom)>=10 and max(hip_rom)<=80
    accepted=accepted and min(coverage)>=.8 and phase_error<.15 and .8<=T<=2.5 and cv<.35
    return dict(accepted=bool(accepted),principal_knee_cycles=list(map(len,peaks)),principal_peak_times_s=[p.tolist()for p in peaks],period_s=T,period_cv=cv,bilateral_phase_error=phase_error,cycle_coverage=min(coverage),knee_ranges_deg=knee_rom.tolist(),hip_ranges_deg=hip_rom.tolist(),autocorrelation=correlations,dominant_spectral_Hz=spectral,peak_counts=secondary)

def self_check():
    t=np.arange(0.,9.,.002);T=2.466
    # A repeatable flexion waveform with a secondary, smaller flexion hump.
    def wave(phase):
        x=(t/T+phase)%1.
        return 65*np.exp(-((x-.25)/.09)**2)+35*np.exp(-((x-.48)/.06)**2)
    hip=np.column_stack((20*np.sin(2*np.pi*t/T),20*np.sin(2*np.pi*t/T+np.pi)))
    anti=measure(t,np.column_stack((wave(0),wave(.5))),hip)
    sync=measure(t,np.column_stack((wave(0),wave(0))),hip)
    fast=measure(t,np.column_stack((35+30*np.sin(2*np.pi*10*t),35-30*np.sin(2*np.pi*10*t))),hip)
    frozen=measure(t,np.ones((len(t),2)),hip)
    result=dict(double_hump_antiphase=anti,synchronous=sync,fast_10Hz=fast,frozen=frozen)
    result['passed']=bool(anti['accepted']and not any(r['accepted']for r in(sync,fast,frozen)))
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    check=self_check()
    if not check['passed']:raise ValueError('Physical-cycle counter self-check failed: '+json.dumps(check))
    z=np.load(a.source/'aligned_trace.npz');t=z['t_state'];keep=(t>=5.)&(t<14.);names=z['joint_names'].tolist();q=z['qpos'][:,z['joint_qpos_addresses']]
    knee=-np.degrees(q[keep][:,[names.index('knee_angle_'+s)for s in('r','l')]])
    hip=np.degrees(q[keep][:,[names.index('hip_flexion_'+s)for s in('r','l')]])
    result=measure(t[keep],knee,hip);result.update(self_check=check,source=str(a.source),scope='Principal sagittal physical-cycle audit only; preserve legacy metric and require separate full duration, knee, frontal, support, and clearance gates')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
