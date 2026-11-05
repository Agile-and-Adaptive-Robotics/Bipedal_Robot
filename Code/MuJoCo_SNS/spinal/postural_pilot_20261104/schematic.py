"""Editable diagram and exact edge ledger of the tested postural subnetwork."""
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse,FancyArrowPatch,Circle,Polygon,FancyBboxPatch
from model import build,NAMES,R

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--fit',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    parameters=json.loads((a.fit/'FIT_COMPLETE.json').read_text())['parameters']if a.fit else {}
    n=build('premotor','inhibitory',parameters=parameters)
    pos=dict(PP_POS=(3.,5.),PP_NEG=(6.,5.),PP_IN_POS=(3.8,3.8),PP_IN_NEG=(5.2,3.8),MN_POS=(3.,1.7),MN_NEG=(6.,1.7),Ia_POS=(0.,4.5),II_POS=(0.,5.3),Ib_POS=(0.,2.8),IbRelay_POS=(1.5,2.8),Ia_NEG=(9.,4.5),II_NEG=(9.,5.3),Ib_NEG=(9.,2.8),IbRelay_NEG=(7.5,2.8))
    labels={name:name.replace('PP_IN','Postural IN').replace('PP_','Posture ').replace('MN_','MN ').replace('IbRelay_','Ib IN ').replace('_',' ')for name in NAMES}
    labels.update(PP_IN_POS='Postural IN\nPOS',PP_IN_NEG='Postural IN\nNEG')
    fig,ax=plt.subplots(figsize=(11,7));ax.set(xlim=(-.9,9.9),ylim=(.15,7.));ax.set_aspect('equal');ax.axis('off');plt.rcParams['svg.fonttype']='none'
    ax.add_patch(FancyBboxPatch((2.2,.8),4.6,4.8,boxstyle='round,pad=.1',fc='#f6f7f8',ec='.5',ls='--',zorder=0))
    patches={}
    for name,(x,y)in pos.items():
        inhibitory=name.startswith(('PP_IN','IbRelay'))
        color='#DDCC77'if inhibitory else'#88CCEE'if name.startswith(('PP_','MN_'))else'#B9DEC7'
        width=1.15 if name.startswith(('PP_IN','IbRelay'))else 1.2
        patch=Ellipse((x,y),width,.58 if name.startswith('PP_IN')else.5,fc=color,ec='.15',lw=1.2,zorder=3)
        patches[name]=patch;ax.add_patch(patch);ax.text(x,y,labels[name],ha='center',va='center',fontsize=10,zorder=4)
    def connect(source,target,inhibitory=False,curvature=0.,style='-'):
        x,y=pos[source];xx,yy=pos[target]
        arrow=FancyArrowPatch((x,y),(xx,yy),patchA=patches[source],patchB=patches[target],connectionstyle='arc3,rad='+str(curvature),arrowstyle='-'if inhibitory else'-|>',mutation_scale=13,shrinkA=2,shrinkB=2,lw=1.2,ec='.18',fc='white',linestyle=style,zorder=1)
        ax.add_patch(arrow)
        if inhibitory:
            # Place a filled circle at the trimmed target end of the actual path.
            path=arrow.get_path();endpoint=path.vertices[-1];ax.add_patch(Circle(endpoint,.055,fc='.12',ec='.12',zorder=2))
    for edge in n.edges:
        s,t=edge['source'],edge['target'];curve=0.
        if s.startswith('PP_IN')and t.startswith('PP_'):curve=.16 if s.endswith('POS')else-.16
        if s.startswith('Ia_')and t.startswith('MN_'):curve=.14 if s.endswith('POS')else-.14
        connect(s,t,edge['reversal_mV']<0,curve,'--'if s.startswith('Ia_')and t.startswith('MN_')else'-')
    for x,side in ((3.,'POS'),(6.,'NEG')):
        ax.add_patch(Polygon([[x-.22,6.3],[x+.22,6.3],[x,5.95]],fc='#228833',ec='.15',zorder=3))
        ax.annotate('',xy=(x,5.3),xytext=(x,5.95),arrowprops=dict(arrowstyle='-|>',ec='.15',fc='white'))
        ax.text(x,6.5,'Descending '+side+' command',ha='center',fontsize=10)
        ax.annotate('Tensile actuator '+side,xy=(x,1.4),xytext=(x,.45),ha='center',fontsize=10,arrowprops=dict(arrowstyle='<-',color='.15'))
    ax.set_title('Tested antagonistic postural premotor subnetwork',fontsize=15,pad=14)
    fig.text(.5,.08,'Open triangle: excitation   |   Filled circle: inhibition   |   Dashed: local Ia to MN',ha='center',fontsize=10)
    fig.text(.5,.045,'POS/NEG indicate measured joint-torque signs. Synaptic encoding span: 0–5 mV relative to rest.',ha='center',fontsize=10)
    for ext in ('svg','pdf','png'):fig.savefig(a.output/('POSTURAL_SUBNETWORK.'+ext),dpi=180,bbox_inches='tight')
    ledger=dict(nodes=[dict(name=name,position=pos[name],tau_s=.01 if name.startswith('MN_')else .02*parameters.get('tau_scale',1.)if name.startswith('PP_')else .02)for name in NAMES],edges=[dict(**edge,e_lo_mV=0.,e_hi_mV=R)for edge in n.edges],parameters=parameters,scope='Actual isolated ankle-bench circuit; no RG or PF populations in this module')
    (a.output/'POSTURAL_SUBNETWORK.json').write_text(json.dumps(ledger,indent=2))
if __name__=='__main__':main()
