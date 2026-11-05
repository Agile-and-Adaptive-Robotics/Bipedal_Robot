"""Measured neural/joint spectral and peak diagnostics from one saved run."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.signal import periodogram,find_peaks

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    z=np.load(a.source/'aligned_trace.npz');t=z['t_state'];mask=(t>=5.)&(t<14.);dt=float(np.median(np.diff(t)))
    rows=[]
    def inspect(name,signal,units):
        y=signal[mask];f,power=periodogram(y,fs=1/dt,window='hann');peaks=find_peaks(y,prominence=max(.05*np.ptp(y),1e-7),distance=max(1,round(.1/dt)))[0]
        spectral=find_peaks(power)[0];spectral=spectral[spectral>0];order=spectral[np.argsort(power[spectral])[::-1]][:4]
        rows.append(dict(channel=name,units=units,range=float(np.ptp(y)),peak_count=len(peaks),median_peak_interval_s=float(np.median(np.diff(t[mask][peaks])))if len(peaks)>2 else None,spectral_peaks_Hz=f[order].tolist(),spectral_power=power[order].tolist()))
    for i,name in enumerate(z['neural_names'].tolist()):
        if name.startswith(('RG_','PF_','MN_')):inspect(name,z['neural_voltage'][:,i],'rest-relative mV')
    for i,name in enumerate(z['joint_names'].tolist()):inspect(name,np.degrees(z['qpos'][:,z['joint_qpos_addresses'][i]]),'deg')
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(dict(source=str(a.source),window_s=[5,14],channels=rows),indent=2))
    print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
