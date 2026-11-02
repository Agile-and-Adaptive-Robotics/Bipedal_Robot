"""Isolated additive, time-bounded pelvis pulse; no historical runner edits."""
import ast,json
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
import numpy as np

class SinglePulse:
    def __init__(self,force_N=0.,onset_s=4.,duration_s=.15,axis='+x'):
        if axis not in ('+x','-x','+y','-y','+z','-z'):raise ValueError('Explicit native-world signed axis required')
        if not np.isfinite(force_N+onset_s+duration_s) or force_N<0 or onset_s<0 or duration_s<=0:raise ValueError('Invalid pulse')
        self.force=float(force_N);self.onset=float(onset_s);self.duration=float(duration_s);self.axis='xyz'.index(axis[-1]);self.sign=-1. if axis[0]=='-' else 1.
    def vector(self,time_s):
        f=np.zeros(3)
        if self.onset<=time_s<self.onset+self.duration:f[self.axis]=self.sign*self.force
        return f
    @contextmanager
    def applied(self,data,body,time_s):
        force=self.vector(time_s);baseline=data.xfrc_applied[body].copy();data.xfrc_applied[body,:3]+=force
        try:yield force
        finally:data.xfrc_applied[body]=baseline

def regression(out):
    runner=Path(__file__).resolve().parent.parent/'runner.py';tree=ast.parse(runner.read_text(encoding='utf-8'))
    node=next(n for n in ast.walk(tree) if isinstance(n,ast.If) and n.lineno>1700 and ast.unparse(n.test)=='push_N > 0.0')
    code=compile(ast.Module(body=[node],type_ignores=[]),str(runner),'exec')
    data=SimpleNamespace(xfrc_applied=np.zeros((1,6)));legacy_data=SimpleNamespace(xfrc_applied=np.zeros((1,6)));fake=SimpleNamespace(mj_name2id=lambda *args:0,mjtObj=SimpleNamespace(mjOBJ_JOINT=0));env=dict(push_N=20.,push_time=4.,push_axis='+x',push_state=[False,None],mujoco=fake,model=SimpleNamespace(jnt_bodyid=[0]),data=legacy_data)
    t=np.arange(4000)*.002;legacy=[];new=[];z=[];pulse=SinglePulse(20.)
    for time_s in t:
        legacy_data.xfrc_applied[0,2]=10.;env['t']=float(time_s);exec(code,env);legacy.append(legacy_data.xfrc_applied[0,0])
        data.xfrc_applied[0]=[0.,0.,10.,1.,2.,3.]
        with pulse.applied(data,0,float(time_s)):new.append(data.xfrc_applied[0].copy())
        z.append(data.xfrc_applied[0].copy())
    new=np.asarray(new);z=np.asarray(z);legacy=np.asarray(legacy);mask=new[:,0]!=0
    transitions=int(np.count_nonzero(np.diff(np.r_[False,mask].astype(int))==1))
    assert transitions==1 and mask.sum()==75 and np.all(new[~mask,0]==0.) and np.all(new[:,2]==10.)
    assert np.array_equal(z,np.tile([0.,0.,10.,1.,2.,3.],(len(t),1)))
    data.xfrc_applied[0]=[7.,8.,9.,1.,2.,3.];baseline=data.xfrc_applied.copy()
    with SinglePulse(0.).applied(data,0,4.05):assert np.array_equal(data.xfrc_applied,baseline)
    assert np.array_equal(data.xfrc_applied,baseline)
    result=dict(runner_source=str(runner),extracted_legacy_if_line=node.lineno,dt_s=.002,legacy_nonzero_steps=int(np.count_nonzero(legacy)),legacy_restarts=int(np.count_nonzero(np.diff(np.r_[False,legacy!=0].astype(int))==1)),new_nonzero_steps=int(mask.sum()),new_pulses=transitions,new_impulse_Ns=float(new[:,0].sum()*.002),other_wrench_preserved=True,zero_default_bit_identical=True,scope='Direct AST execution of existing runner push block versus isolated additive pulse; no production mutation.')
    out.write_text(json.dumps(result,indent=2));print(json.dumps(result))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();regression(a.out)
