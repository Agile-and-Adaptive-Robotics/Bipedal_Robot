"""Small layer-specific current perturbations of the verified Hunt air model."""
import argparse,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--amplitude',type=float,default=.5);p.add_argument('--start',type=float,default=7.1);p.add_argument('--duration',type=float,default=.1);p.add_argument('--sham',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);rows=[]
    if not 0<a.amplitude<=2.:raise ValueError('Bounded current amplitude required')
    cases=[('RG','RG_F_r',a.amplitude),('PF','PF_HIP-F_r',a.amplitude),('MN','MN_iliacus_r',a.amplitude),('PF_inhibition','PF_HIP-F_r',-a.amplitude)]
    if a.sham:cases.insert(0,('Baseline','RG_F_r',0.))
    for label,cell,amplitude in cases:
        folder=a.output/label
        with(a.output/(label+'.log')).open('a',encoding='utf-8')as log:
            run=subprocess.run([sys.executable,str(HERE/'body_probe.py'),'--reflex-fit',str(a.fit),'--stage','2','--expression','.1','--stim-cell',cell,'--stim-nA',str(amplitude),'--stim-start',str(a.start),'--stim-duration',str(a.duration),'--output',str(folder)],stdout=log,stderr=subprocess.STDOUT,timeout=600,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        path=folder/'body.result.json';row=dict(layer=label,cell=cell,current_nA=amplitude,start_s=a.start,duration_s=a.duration,source=str(folder),exit_code=run.returncode)
        if path.exists():
            result=json.loads(path.read_text());row.update(complete=result['assessment']['complete'],knee_safe=result['assessment']['knee_extension_gate'],assessment=result['assessment'])
        elif(folder/'STOP.json').exists():row['stop']=json.loads((folder/'STOP.json').read_text())
        rows.append(row);(a.output/'STATUS.json').write_text(json.dumps(rows,indent=2))
    (a.output/'COMPLETE.json').write_text(json.dumps(rows,indent=2))
if __name__=='__main__':main()
