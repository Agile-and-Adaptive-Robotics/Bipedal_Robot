"""Build the Appendix D ground-trial figure from paired, matched reruns."""
from pathlib import Path
import argparse
import hashlib
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

def load_case(folder):
    with np.load(folder/"output.npz",allow_pickle=False) as z,np.load(folder/"true_activation.npz",allow_pickle=False) as activation:
        nn=z["neuro_names"].tolist();jn=z["key_joints"].tolist();expected=["vas_lat_r","semimem_r"]
        assert activation["muscle_names"].tolist()==expected
        label=z["t"];dt=float(np.median(np.diff(label)))
        state=activation["state_time"];offset=len(label)-len(state)
        np.testing.assert_allclose(state,label[offset:]+dt,atol=1e-9)
        use=(state>=10.-1e-9)&(state<=12.5+1e-9)
        result={"time":state[use]}
        for field,channels in [("rg",["RG_E_r","RG_F_r"]),("pf",["PF_E1_r","PF_F1_r"])]:
            result[field]=z["neuro"][offset:][use][:,[nn.index(v) for v in channels]]
        result["activation"]=activation["activation"][use]
        for label,j in [("hip","hip_flexion"),("knee","knee_angle"),("ankle","ankle_angle")]:
            result[label]=z["q"][offset:][use][:,[jn.index(j+"_r"),jn.index(j+"_l")]]
        return result

def render(package):
    with np.load(package,allow_pickle=False) as data:
        plt.rcParams.update({"font.family":"Arial","font.size":12.5,"axes.labelweight":"bold","axes.linewidth":1.15,"xtick.major.width":1.15,"ytick.major.width":1.15,"pdf.fonttype":42,"svg.fonttype":"none","path.simplify":False})
        fig,ax=plt.subplots(6,2,figsize=(7.2,8.5),sharex=True,sharey="row")
        fig.subplots_adjust(left=.15,right=.985,top=.867,bottom=.073,hspace=.20,wspace=.19)
        styles=[("#D55E00","-"),("#0072B2","--")]
        for col,title in enumerate(["Nonspiking","Spiking (untuned)"]):
            ax[0,col].set_title(title,fontweight="bold",pad=7)
            t=data[f"time_{col}"]
            for row,(field,ylab) in enumerate([("rg","RG (mV)"),("pf","PF (mV)"),("activation","Muscle\nactivation"),("hip","Hip (deg)"),("knee","Knee (deg)"),("ankle","Ankle (deg)")]):
                a=ax[row,col];values=data[f"{field}_{col}"]
                for k,(color,style) in enumerate(styles):
                    a.plot(t,values[:,k],color=color,ls=style,lw=1.45)
                a.spines["top"].set_visible(False);a.spines["right"].set_visible(False)
                a.yaxis.set_major_locator(MaxNLocator(nbins=3))
                a.tick_params(axis="both",labelsize=12.5)
                a.set_xlim(10.,12.5)
                if col==0:a.set_ylabel(ylab)
                if row==5:a.set_xlabel("Time (s)");a.set_xticks([10.,11.,12.])
        legend=fig.legend([Line2D([],[],color=c,ls=ls,lw=1.7) for c,ls in styles],["Extensor / Right","Flexor / Left"],loc="upper center",bbox_to_anchor=(.55,.995),ncol=2,frameon=True,fancybox=False,title="Legend",handlelength=2.1)
        legend.get_frame().set_edgecolor("black");legend.get_frame().set_linewidth(1.0);legend.get_title().set_fontweight("bold")
        path=package.parents[2]/"spiking_ground_neuromechanics_short.pdf"
        fig.savefig(path,metadata={"Title":"Nonspiking and spiking Gait2392 ground trials","CreationDate":None,"ModDate":None})
        plt.close(fig)
        return path

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--runs",type=Path,help="Folder with independent figpair_s3k and figpair_spiking runs.")
    parser.add_argument("--package",type=Path,required=True)
    a=parser.parse_args()
    if a.runs:
        a.package.parent.mkdir(parents=True,exist_ok=True)
        arrays={};meta={}
        for k,name in enumerate(["s3k","spiking"]):
            d=a.runs/("figpair_"+name);out=load_case(d)
            arrays.update({key+"_"+str(k):value for key,value in out.items()})
            meta[name]={"sources":{str(d/f):hashlib.sha256((d/f).read_bytes()).hexdigest() for f in ["output.npz","true_activation.npz","result.json"]},"sample_count":len(out["time"]),"start_s":float(out["time"][0]),"end_s":float(out["time"][-1]),"channels":{"rg":["RG_E_r","RG_F_r"],"pf":["PF_E1_r","PF_F1_r"],"activation":["vas_lat_r","semimem_r"]}}
        np.savez_compressed(a.package,**arrays)
        (a.package.parent/"provenance.json").write_text(json.dumps(meta,indent=2))
    print(render(a.package))
if __name__=="__main__":main()
