"""Short drive-response mapping of the frozen rhythm core, without body physics."""
import argparse,json,os,sys,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
from replay_adapter import Replay

def main():
    p=argparse.ArgumentParser();p.add_argument('--timing',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    cfg=json.loads((a.timing/'MANIFEST.json').read_text());fixed=json.loads((a.timing/'FINAL.json').read_text())['parameters'];r=Replay(cfg);r.u.fill(0.);rows=[]
    fixed.update(Ia_RG=0.,II_RG=0.,Ib_RG=0.,Ia_PF=0.,II_PF=0.,Ib_PF=0.)
    cases=[(v,v)for v in(0.,.15,.35,.55,.75,1.)]+[(.75,.25),(1.,.25),(1.,.5)]
    for flex,ext in cases:
        r.manifest.update(command_F=flex,command_E=ext);result=r.evaluate(fixed,time.time()+180.)
        rows.append(dict(F=flex,E=ext,**result));(a.output/'RESULTS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(dict(F=flex,E=ext,**result)),flush=True)
    (a.output/'MANIFEST.json').write_text(json.dumps(dict(source_timing=str(a.timing),afferents='zero inputs and zero central gains',configuration=cfg,scope='Neural rhythm-core operating map; not body or walking acceptance'),indent=2))
if __name__=='__main__':main()
