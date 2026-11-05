"""Scientific traces and replay animation from saved postural-bench states."""
import argparse,json,os,sys
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import mujoco
from model import XML,NAMES

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none'})
    fig,axes=plt.subplots(4,1,figsize=(8.3,7),sharex=True,constrained_layout=True)
    colors=('#4477AA','#228833','#CC6677')
    for policy,color in zip(('off','inhibitory','support_positive'),colors):
        z=np.load(a.input/('premotor_'+policy+'.npz'));log=z['mechanics'];neural=z['neural'];keep=(log[:,0]>=1.4)&(log[:,0]<=4.3);t=log[keep,0]
        axes[0].plot(t,np.degrees(log[keep,1]),label=policy.replace('_',' '),color=color)
        if policy=='inhibitory':
            for name,style in (('PP_POS','-'),('PP_NEG','--')):axes[1].plot(t,neural[keep,NAMES.index(name)],style,label=name,color='#4477AA' if style=='-'else'#EE6677')
            for name,style in (('MN_POS','-'),('MN_NEG','--')):axes[2].plot(t,neural[keep,NAMES.index(name)],style,label=name,color='#4477AA' if style=='-'else'#EE6677')
            axes[3].plot(t,log[keep,5],label='positive torque actuator',color='#4477AA');axes[3].plot(t,log[keep,6],'--',label='negative torque actuator',color='#EE6677')
    for ax in axes:
        ax.axvspan(2.,2.05,color='#DDCC77',alpha=.5);ax.grid(alpha=.2);ax.legend(loc='upper right',fontsize=8)
    axes[0].set_ylabel('Ankle angle (deg)');axes[1].set_ylabel('Premotor (mV)');axes[2].set_ylabel('MN (mV)');axes[3].set_ylabel('Tension (N)');axes[3].set_xlabel('Time (s)')
    fig.suptitle('Isolated antagonistic postural subnetwork: 8 N m push for 50 ms')
    for ext in ('png','pdf','svg'):fig.savefig(a.output/('POSTURAL_RESPONSE.'+ext),dpi=180)
    plt.close(fig)
    z=np.load(a.input/'premotor_inhibitory.npz');log=z['mechanics']
    m=mujoco.MjModel.from_xml_string(XML);d=mujoco.MjData(m);renderer=mujoco.Renderer(m,height=480,width=640)
    try:
        # 30 fps, replaying actual sampled physical states; no second dynamics run.
        frames=[]
        for t in np.arange(0.,5.,1/30):
            d.qpos[0]=np.interp(t,log[:,0],log[:,1]);mujoco.mj_forward(m,d);renderer.update_scene(d,camera='side');frames.append(Image.fromarray(renderer.render()))
        frames[0].save(a.output/'POSTURAL_ANKLE_BENCH.gif',save_all=True,append_images=frames[1:],duration=33,loop=0)
        frames[0].save(a.output/'POSTURAL_ANKLE_BENCH.png')
    finally:
        if hasattr(renderer,'close'):renderer.close()
    (a.output/'README.md').write_text('Isolated 10 kg ankle bench; fixed base, no harness. The yellow plot band is an 8 N m push for 50 ms. The neural subnetwork receives Ia/II/Ib inputs, with an explicit inhibitory Ib alternative. The GIF replays saved states and does not show Gait2392 or robot walking.\n')
if __name__=='__main__':main()
