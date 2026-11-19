"""Rebuild the chat's bilateral syn6 overlay with publication-size vector labels."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.signal import find_peaks

here = Path(__file__).resolve().parent
diss = here.parents[1]
repo = diss.parents[2]
spinal = repo/"Code/MuJoCo_SNS/spinal"
tools = spinal/"controller_profiles_20261104"
sys.path.insert(0,str(tools))
from kinematic_response import reference, JOINTS, GRID

p=argparse.ArgumentParser(description=__doc__)
p.add_argument("--records",type=Path,default=Path("D:/Users/Ben Bolen/ControllerTraining_20261104"))
p.add_argument("--use-saved",action="store_true",help="Render from the portable, verified plot data.")
args=p.parse_args()
figdir=diss/"Figures/60-sns"
component=figdir/"Components/syn6_stage2_overlay"
component.mkdir(parents=True,exist_ok=True)
arrays=component/"plot_data.npz"
if not args.use_saved:
    selected=args.records/"syn6_s3_single_20261105/left_right_figures"
    record=selected/"OVERLAY_METRICS.json"
    expected=json.loads(record.read_text())
    original_img=selected/"KINEMATIC_OVERLAY.png"
    assert hashlib.sha256(original_img.read_bytes()).hexdigest()=="f66123d0164bf021c966c4edc673d455bf6e104ecb8622cf9adab19f3ee19a5d"
    run=args.records/"syn6_stage2_regularize_20261105/s2_latch"
    assert Path(expected["run"]).name==run.name
    trace=run/"capture/aligned_trace.npz"
    with np.load(trace,allow_pickle=False) as z:
        t=z["t_state"]
        jnames=z["joint_names"].tolist()
        idx=[z["joint_qpos_addresses"][jnames.index(name)] for name in JOINTS]
        angles=np.degrees(z["qpos"][:,idx])
        foot=z["foot_force_anatomical_N"][:,:,1]
    assert np.max(np.abs(foot))==0
    start=max(7.,json.loads((run/"CADENCE.json").read_text())["full_native_since_neural_s"])
    candidates=find_peaks(-angles[:,1],distance=round(.48/np.median(np.diff(t))),prominence=12.)[0]
    edges=t[candidates[(t[candidates]>=start)&(t[candidates]<=14.)]]
    cycles=np.array([np.column_stack([np.interp(a+GRID*(b-a),t,angles[:,i]) for i in range(6)])
                     for a,b in zip(edges[:-1],edges[1:])])
    assert len(cycles)==7
    np.testing.assert_allclose(edges[[0,-1]],expected["cycle_window_s"])
    ong=spinal/"gait_refs/ong_speed_125.npz"
    refs=[reference(),reference(ong)]
    observed_ranges=np.ptp(cycles.mean(axis=0),axis=0)
    for ref,saved in zip(refs,expected["benchmarks"]):
        np.testing.assert_allclose(observed_ranges,[j["simulation_ROM_deg"] for j in saved["joints"]])
        np.testing.assert_allclose(np.ptp(ref["cycle_deg"],axis=0),[j["reference_ROM_deg"] for j in saved["joints"]])
    bar_starts=[[( (r["reference_left_loading_phase"] if s else 0.)-r["right_knee_anchor_fraction"])%1. for s in range(2)] for r in refs]
    np.savez_compressed(arrays,phase=GRID,cycles=cycles,ref=np.array([r["cycle_deg"] for r in refs]),
                        bar_starts=bar_starts,duties=[r["contact_duty"] for r in refs],edges=edges)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    sources=[trace,record,original_img,run/"body.result.json",run/"CADENCE.json",
             tools/"stage2_timing_figures.py",tools/"kinematic_response.py",ong]
    src={str(path):sha(path) for path in sources}
    for ref in refs:
        src.update(ref["source_sha256"])
    audit=dict(experiment="Completed afferented, suspended Gait2392 Stage 2 with an S2 shutoff latch",
               original_figure=str(original_img),sources_sha256=src,cycles=len(cycles),window_s=expected["cycle_window_s"],
               aligning_event="Measured maximum right-knee flexion shared by all six channels. No per-joint alignment.",
               foot_vertical_load_N_peak=float(np.abs(foot).max()),PF_bursts=expected["pf_onset_counts_by_rg_cycle"],
               ranges_of_mean_LL_RL_deg={"syn6":observed_ranges[[4,1]].tolist(),
                    "Subject01":np.ptp(refs[0]["cycle_deg"],axis=0)[[4,1]].tolist(),
                    "Ong":np.ptp(refs[1]["cycle_deg"],axis=0)[[4,1]].tolist()},
               mean_of_individual_excursions_left_right_deg=np.ptp(cycles[:,:,[4,1]],axis=1).mean(axis=0).tolist(),
               plot_data_sha256=sha(arrays),phase_samples=len(GRID))
    np.testing.assert_allclose(audit["mean_of_individual_excursions_left_right_deg"],[48.6,49.6],atol=.055)
    (here/"SYN6_OVERLAY.json").write_text(json.dumps(audit,indent=2)+"\n",encoding="utf-8")

with np.load(arrays) as z:
    x=100*z["phase"]
    cycles=z["cycles"]
    refs=z["ref"]
    starts=z["bar_starts"]
    duties=z["duties"]
plt.rcParams.update({"font.family":"Arial","font.size":10.5,"xtick.labelsize":10,"ytick.labelsize":10,
                     "axes.labelsize":10.5,"axes.titlesize":11,"pdf.fonttype":42,"path.simplify":False})
fig,ax=plt.subplots(4,2,figsize=(6.35,5.8),sharex=True,sharey="row",
                    gridspec_kw={"height_ratios":[1,1,1,.46]})
fig.subplots_adjust(left=.148,right=.962,top=.874,bottom=.112,hspace=.17,wspace=.095)
blue,gray,orange="#2266A5","#262626","#B44B13"
for column,(side,leg) in enumerate([("Left",1),("Right",0)]):
    ax[0,column].set_title(side,weight="bold",pad=3)
    for row,joint in enumerate(("Hip","Knee","Ankle")):
        panel=ax[row,column]
        i=3*leg+row
        for series in cycles[:,:,i]:
            panel.plot(x,series,color=blue,alpha=.24,lw=.62)
        panel.plot(x,cycles[:,:,i].mean(axis=0),color=blue,lw=1.5)
        panel.plot(x,refs[0,:,i],color=gray,ls="--",lw=1.32)
        panel.plot(x,refs[1,:,i],color=orange,ls=":",lw=1.7)
        if column==0:
            panel.set_ylabel(joint+" (°)")
        panel.spines[["right","top"]].set_visible(False)
        panel.tick_params(length=2.6)
    panel=ax[3,column]
    for j,color in enumerate([gray,orange]):
        y=1-j
        onset=starts[j,leg]
        duty=duties[j,leg]
        panel.broken_barh([(0,100)],(y-.24,.48),facecolors="#e0e0e0")
        parts=[(100*onset,100*min(duty,1-onset))]
        if onset+duty>1:
            parts.append((0,100*(onset+duty-1)))
        panel.broken_barh(parts,(y-.24,.48),facecolors=color)
    panel.set(yticks=[1,0],yticklabels=["Subject01","Ong"],ylim=(-.48,1.48),xlim=(0,100),xticks=[0,25,50,75,100])
    panel.spines[["right","top","left"]].set_visible(False)
    panel.tick_params(axis="y",length=0)
    panel.tick_params(axis="x",length=2.6)
fig.legend(handles=[Line2D([],[],color=blue,lw=1.7,label="Suspended syn6 mean"),
                    Line2D([],[],color=gray,ls="--",lw=1.5,label="Subject01 walking"),
                    Line2D([],[],color=orange,ls=":",lw=1.8,label="Ong simulated walking")],
           loc="upper center",bbox_to_anchor=(.5,1.005),ncol=2,frameon=False,fontsize=10.2,
           columnspacing=1.4,handlelength=2.1)
fig.text(.565,.018,"Cycle from maximum right-knee flexion (%)",ha="center",fontsize=10.5)
fig.canvas.draw()
for panel in ax.flat:
    for tick in panel.get_xticklabels()+panel.get_yticklabels():
        if tick.get_visible() and tick.get_text():
            assert fig.bbox.contains(*tick.get_window_extent().get_points()[0])
            assert fig.bbox.contains(*tick.get_window_extent().get_points()[1])
fig.savefig(figdir/"syn6_suspended_kinematics.pdf",metadata={"Title":"Bilateral movements of the suspended syn6 model and two walking references","Creator":"syn6_overlay.py"})
plt.close(fig)
(figdir/"syn6_suspended_kinematics.alt.txt").write_text("Left and right hips, knees, and ankles of suspended syn6 are compared with recorded Subject01 walking and walking simulated by Ong and colleagues. Light curves are the seven individual suspended cycles. The blue mean knee excursions are 19 degrees on the left and 34 degrees on the right, against about 70 to 72 degrees in the reference curves. All joints share the time origin of maximum right-knee flexion. Beneath each leg, dark bars mark reference stance and pale gray bars mark swing. Neither suspended foot carried load. Knee flexion is negative.\n",encoding="utf-8")
