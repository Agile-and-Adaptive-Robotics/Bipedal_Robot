"""Compact display-only revision using exact saved arrays and shared BPA geometry.

Derived from the prior figure renderer; no experiments, interpolation, smoothing,
phase shifting, or decimation. Recreate TIMING_trace_arrays.npz with TIMING_audit.py.
"""
from pathlib import Path
import json
import os
import hashlib

ROOT = next(p for p in Path(__file__).resolve().parents if p.name == "Bipedal_Robot")
HERE = Path(__file__).parent
os.environ["MPLCONFIGDIR"] = str(ROOT / "tmp/ch6_authoritative_20261116/mplconfig")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle, Rectangle
import numpy as np

TRACE=ROOT/'tmp/ch6_authoritative_20261116/TIMING_trace_arrays.npz'
Z=np.load(TRACE,allow_pickle=False)
keep=(Z['knee_t']>=29)&(Z['knee_t']<=40)
DATA={'knee':{'t':Z['knee_t'][keep],'P':Z['knee_P'][keep],'C':Z['knee_C'][keep]},'cup':[]}
for condition in ['on','off']:
    run={}
    for key in ['th','mCup','F','A','Ft','At']:
        stem='cup_graded_'+condition+'_'+key
        run[key]={'t':Z[stem+'_t'],'y':Z[stem+'_y']}
    DATA['cup'].append(run)
# Exact Python expression of the shared sns_icon_geom.m BPA branch.
xb=np.linspace(-.70,.70,49);rb=.185+.115*(1-(xb/.70)**2)
def radius(x):return .185+.115*max(1-(x/.70)**2,0)
braid=[]
for i in range(7):
    x1=-.64+.175*i;x2=min(x1+.30,.66)
    braid.extend([[x1,-.90*radius(x1),x2,.90*radius(x2)],[x1,.90*radius(x1),x2,-.90*radius(x2)]])
DATA['bpa']={'bodyX':np.r_[xb,xb[::-1]],'bodyY':np.r_[rb,-rb[::-1]],'braid':braid,
    'fill':[.94,.76,.74],'cap':[.55,.55,.55],'capL':[[-.70,-.20],[-.90,-.20],[-.90,.20],[-.70,.20]],
    'capR':[[.70,-.20],[.90,-.20],[.90,.20],[.70,.20]]}
OUT = HERE/'CAD_sns_compact'
OUT.mkdir(exist_ok=True)
BLUE = "#4477AA"
RED = "#EE6677"
plt.rcParams.update({"font.family": "Arial", "font.size": 10.5, "axes.labelsize": 10.5,
                     "axes.titlesize": 11, "axes.linewidth": 1,
                     "pdf.fonttype": 42, "ps.fonttype": 42,
                     "path.simplify": False, "savefig.facecolor": "white"})

def figure(w, h):
    return plt.figure(figsize=(w / 2.54, h / 2.54), facecolor="white")

def save(fig, stem, alt):
    fig.savefig(OUT / (stem + ".pdf"))
    fig.savefig(OUT / (stem + ".png"), dpi=200)
    (OUT / (stem + ".alt.txt")).write_text(alt + "\n", encoding="utf-8")
    plt.close(fig)
    print("Wrote", stem)

def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")
    ax.grid(True, color="0.88", lw=.55)
    ax.set_axisbelow(True)

def panel(ax,letter):
    ax.text(-.185,1.03,letter,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom')

def glyph(ax, items):
    for p in items:
        if p["kind"] == "patch":
            ax.add_patch(Polygon(np.column_stack([p["x"], p["y"]]),
                                 facecolor=p["colour"], edgecolor="none"))
        elif p["kind"] == "line":
            ax.plot(p["x"], p["y"], color=p["colour"], lw=.9)
        else:
            ax.text(0, 0, p["text"], ha="center", va="center", weight="bold")
    ax.set_aspect("equal")
    ax.set_axis_off()

def library():
    names = ["NonSpikingNeuron", "NonSpikingSynapse", "SynSum", "SpikingLIFNeuron",
             "IaMuscleSpindle", "IbGolgiTendon", "MuscleActivation", "BioMuscle", "BPAForce",
             "BPA_10mm", "BPA_20mm", "BPA_40mm"]
    labels = ["Nonspiking neuron", "Graded synapse", "Synaptic summation", "Spiking LIF neuron",
              "Ia spindle encoder", "Ib force encoder", "Activation dynamics", "Biological muscle",
              "Simple force element", "10 mm BPA", "20 mm BPA", "40 mm BPA"]
    subs = ["Membrane voltage\nand drive", "Excitation or\ninhibition", "Conductance pairs",
            "Threshold and reset", "Length and\nstretch velocity", "Tendon-force signal", "First-order filter",
            "Thelen-style\nformulation", "Reduced-order\nforce law", "Pressure and length", "Pressure and length", "Pressure and length"]
    icons = {p["name"]: p["items"] for p in DATA["icons"]}
    fig = figure(14, 15.2)
    for k, name in enumerate(names):
        row, col = divmod(k, 3)
        cx, top = .025 + col * .327, .99 - row * .247
        fig.add_artist(Rectangle((cx, top - .235), .309, .23, transform=fig.transFigure,
                                 edgecolor=".8", facecolor="white", lw=.6, zorder=-1))
        muscle = "BPA" in name or name == "BioMuscle"
        if muscle:
            ax = fig.add_axes([cx + .037, top - .108, .235, .079])
            ax.set_xlim(-1.05, 1.05); ax.set_ylim(-.4, .4)
        else:
            ax = fig.add_axes([cx + .071, top - .125, .165, .116])
            ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1)
        glyph(ax, icons[name])
        fig.text(cx + .1545, top - .155, labels[k], ha="center", va="center", weight="bold", fontsize=10.5)
        fig.text(cx + .1545, top - .201, subs[k], ha="center", va="center", fontsize=10, linespacing=1.12)
    save(fig, "sns_library_20261116", "Twelve current SNS-Simscape block glyphs matching the selection in defense slide 69, rendered from the actual mask drawing primitives and corrected shared side-view actuator geometry. Graded and spiking neurons have different voltage waveforms; sensory capsules, activation pentagon, synapse summation, pennate biological muscle, simple force element, and 10, 20, and 40 mm BPA models are identified.")

def edge(ax, pts, exc, color="black", lw=1.2):
    pts = np.array(pts, float)
    ax.plot(pts[:, 0], pts[:, 1], "-" if exc else "--", color=color, lw=lw, zorder=1)
    u = pts[-1] - pts[-2]; u /= np.linalg.norm(u)
    n = np.array([-u[1], u[0]])
    q = pts[-1]
    if exc:
        ax.add_patch(Polygon([q - .17*u, q - .12*n, q + .12*n],
                             facecolor="white", edgecolor="black", lw=1, zorder=3))
    else:
        ax.add_patch(Circle(q, .095, facecolor="black", edgecolor="black", zorder=3))

def neuron(ax, p, r, label, color):
    ax.add_patch(Circle(p, r, facecolor=color, edgecolor="black", lw=1.3, zorder=4))
    ax.text(*p, label, ha="center", va="center", weight="bold", zorder=5)

def circuit(ax, names, heading):
    ax.set_xlim(0, 14); ax.set_ylim(0, 5.9); ax.set_axis_off(); ax.set_aspect("equal")
    ax.text(.1, 5.65, heading, weight="bold", fontsize=11)
    sx, sy, my, mx = [2.0, .8, 10.0, 12.8], 4.5, 1.45, [3, 11]
    edge(ax, [(sx[0], sy-.38), (mx[0]-.22, my+.45)], True, BLUE)
    edge(ax, [(sx[1]+.38, 2.8), (mx[0]-.46, my)], False, RED)
    edge(ax, [(sx[2], sy-.38), (mx[1]-.22, my+.45)], True, BLUE)
    edge(ax, [(sx[3]-.38, 2.8), (mx[1]+.46, my)], False, RED)
    edge(ax, [(sx[0]+.38, sy), (6.6, sy), (6.6, my), (mx[1]-.46, my)], False, BLUE)
    edge(ax, [(sx[2]-.38, sy), (7.5, sy), (7.5, my-.8), (mx[0], my-.8), (mx[0], my-.46)], False, BLUE)
    for k, x in enumerate(sx):
        y = sy if k in (0, 2) else 2.8
        neuron(ax, (x, y), .38, "Ia" if k in (0, 2) else "Ib", "#F0F7F0")
        ax.text(x, y+.6, names[k//2], ha="center")
    for k, x in enumerate(mx):
        neuron(ax, (x, my), .46, "MN", "white")
        ax.text(x+(-2.0 if k == 0 else 2.0), my-.55, names[k], ha="center")
    edge(ax, [(.3, .25), (1, .25)], True)
    ax.text(1.25, .25, "excitation", va="center", fontsize=10)
    edge(ax, [(4.5, .25), (5.2, .25)], False)
    ax.text(5.45, .25, "inhibition", va="center", fontsize=10)
    ax.text(9.4, .25, "Ia: stretch   Ib: force", va="center", fontsize=10)

def knee():
    s = DATA["knee"]; t = np.array(s["t"]); p = np.array(s["P"]); c = np.array(s["C"])
    fig = figure(14, 15)
    circuit(fig.add_axes([.025, .645, .95, .335]), ["Extensor", "Flexor"], "A  Implemented antagonist feedback")
    ax = fig.add_axes([.16, .51, .78, .105]); sel = (t >= 29.85) & (t <= 30.25)
    ax.plot(t[sel]-30, p[sel, 8], color=".2", lw=1.5); ax.set(xlim=(-.15,.25), ylim=(-.1,2.3))
    ax.set_ylabel("Ia current\n(nA)"); style(ax)
    panel(ax,'B')
    x = t-30
    ax = fig.add_axes([.16, .295, .78, .14])
    ax.plot(x, -np.rad2deg(p[:, 0]), color=BLUE, lw=1.5, label="Stimulus")
    ax.plot(x, -np.rad2deg(c[:, 0]), "--", color=".3", lw=1.5, label="No stimulus")
    ax.axvline(0, color=".5", ls=":", lw=1); ax.set_xlim(-1,10)
    ax.set_ylabel("Knee angle\n(deg)"); panel(ax,'C')
    ax.legend(loc="lower right", frameon=False, fontsize=10); style(ax)
    ax.text(.99,.94, "Peak +1.09 deg at 111 ms", transform=ax.transAxes, ha="right", va="top", fontsize=10)
    ax = fig.add_axes([.16, .08, .78, .125])
    ax.plot(x, p[:, 2], color=BLUE, lw=1.5, label="Extensor")
    ax.plot(x, p[:, 3], "--", color=RED, lw=1.5, label="Flexor")
    ax.axvline(0, color=".5", ls=":", lw=1); ax.set_xlim(-1,10)
    ax.set_ylabel("MN voltage\n(mV)"); ax.set_xlabel("Time from pulse onset (s)")
    panel(ax,'D')
    ax.legend(loc="upper right", frameon=False, fontsize=10); style(ax)
    save(fig, "sns_knee_reflex_authoritative", "A reduced pathway diagram with all six sensory-to-motor synapses accompanies the saved 2 nA, 30 ms extensor Ia current pulse, knee angle with a matched no-stimulus control, and motor-neuron voltages. The knee extends 1.09 degrees, peaks at 111 ms, and recovers. These are the verified graded-ramp simple-force recordings, not a physiological tendon-impact model or the calibrated Festo force block.")

def bpa(ax, p1, p2, width):
    g=DATA["bpa"]; p1=np.array(p1); p2=np.array(p2); u=(p2-p1)/np.linalg.norm(p2-p1); n=np.array([-u[1],u[0]]); center=(p1+p2)/2
    def trans(x,y):
        return center + np.array(x)[:,None]/1.8*np.linalg.norm(p2-p1)*u + np.array(y)[:,None]*width/.6*n
    xy=trans(g["bodyX"],g["bodyY"]); ax.add_patch(Polygon(xy,facecolor=g["fill"],edgecolor="black",lw=.65,zorder=5))
    for row in g["braid"]:
        xy=trans([row[0],row[2]],[row[1],row[3]]); ax.plot(xy[:,0],xy[:,1],color="black",lw=.45,zorder=6)
    for key in ["capL","capR"]:
        cap=np.array(g[key]); xy=trans(cap[:,0],cap[:,1]); ax.add_patch(Polygon(xy,facecolor=g["cap"],edgecolor="black",lw=.5,zorder=7))

def arm(ax,on,off):
    ax.set(xlim=(-.17,.44),ylim=(-.44,.05)); ax.set_aspect("equal"); ax.set_axis_off()
    sh=np.array([-.06,0]); el=np.array([-.06,-.27])
    ax.plot([sh[0],el[0]],[sh[1],el[1]],color=".6",lw=12,solid_capstyle="round")
    hand=el+.3*np.array([np.cos(np.deg2rad(on)),-np.sin(np.deg2rad(on))])
    offhand=el+.3*np.array([np.cos(np.deg2rad(off)),-np.sin(np.deg2rad(off))])
    ax.plot([el[0],hand[0]],[el[1],hand[1]],color=BLUE,lw=8,solid_capstyle="round")
    ax.plot([el[0],offhand[0]],[el[1],offhand[1]],"--",color=RED,lw=1.5)
    ax.plot(*el,"o",mfc="white",mec="black",ms=7,zorder=4)
    bpa(ax,[-.025,-.045],[.005,-.255],.028); bpa(ax,[-.095,-.045],[-.09,-.28],.022)
    ax.add_patch(Polygon(hand+[[-.035,0],[.035,0],[.04,.1],[-.04,.1]],facecolor="#EEF4FF",edgecolor=".3",lw=1))
    ax.add_patch(Polygon(hand+[[-.035,.005],[.035,.005],[.038,.075],[-.038,.075]],facecolor="#A6D9F2",edgecolor="none"))
    ax.text(.055,-.07,"20 mm BPAs",fontsize=10)
    ax.text(.09,-.405,f"ON {on:.2f} deg\nOFF {off:.2f} deg",ha="center",fontsize=10)

def cup():
    r=DATA["cup"]
    def a(k,key): return np.array(r[k][key]["t"]), np.array(r[k][key]["y"])
    on=np.rad2deg(a(0,"th")[1][-1]); off=np.rad2deg(a(1,"th")[1][-1])
    fig=figure(14,15)
    circuit(fig.add_axes([.04,.724,.92,.258]),["Biceps","Triceps"],"A  Implemented loaded-arm feedback")
    ax=fig.add_axes([.055,.47,.405,.21]); arm(ax,on,off); ax.set_title("B  Arm geometry",loc='left',pad=4)
    ax=fig.add_axes([.61,.535,.32,.15]); ax.plot(*a(0,"mCup"),color=".2",lw=1.5)
    ax.set(xlim=(0,10),ylim=(0,.66),ylabel="Added mass (kg)"); panel(ax,'C'); style(ax)
    ax=fig.add_axes([.17,.285,.77,.155])
    for k,col,ls,lab in [(0,BLUE,"-","Feedback on"),(1,RED,"--","Feedback off")]:
        t,y=a(k,"th"); ax.plot(t,-np.rad2deg(y),ls,color=col,lw=1.5,label=lab)
    ax.set(xlim=(0,10),ylim=(-9,.5),ylabel="Forearm elevation\n(deg)")
    panel(ax,'D'); ax.legend(loc="lower left",frameon=False,fontsize=10); style(ax)
    ax=fig.add_axes([.17,.075,.77,.13]); ax.plot(*a(0,"A"),color=BLUE,lw=1.5,label="Biceps"); ax.plot(*a(0,"At"),"--",color=RED,lw=1.5,label="Triceps")
    ax.set(xlim=(0,10),ylim=(0,.55),ylabel="Pressure command\n(normalized)",xlabel="Time (s)")
    panel(ax,'E'); ax.legend(loc="lower left",frameon=False,fontsize=10); style(ax)
    save(fig,"sns_cup_reflex_authoritative","All six actual feedback connections are shown, with separate biceps and triceps Ib autogenic inhibition. Added mass rises at 0.0625 kg/s throughout the recorded 10-second trial, ending at 0.625 kg. The feedback-on arm sags 2.00 degrees compared with 7.85 degrees for feedback off. Side-view actuator glyphs use current shared vertices. The last panel shows normalized pressure commands for the antagonist Festo 20 mm force blocks; valve/air-flow dynamics are not modeled.")
    return {"cup_final_on_deg":on,"cup_final_off_deg":off,"cup_final_load_kg":a(0,"mCup")[1][-1]}

if __name__ == "__main__":
    knee(); metrics=cup()
    source=[TRACE,ROOT/'Code/Matlab/SNS_Simscape/sns_icon_geom.m',ROOT/'Code/Matlab/SNS_Simscape/demos/sns_build_demo.m',ROOT/'Code/Matlab/SNS_Simscape/demos/sns_build_beer_demo.m',ROOT/'Code/Matlab/SNS_Simscape/results/knee_transient_reflex_20261028/replay_summary.mat',ROOT/'Code/Matlab/SNS_Simscape/results/sns_beer_results.mat',HERE/'TIMING_audit.py']
    old_renderer=HERE.parent/'ch6_revision_20261116/figscripts/render_sns_graphics.py'
    metrics.update(renderer='matplotlib Agg',figure_size_cm=[14,15],minimum_font_pt=10,
       prior_renderer=str(old_renderer),prior_renderer_sha256=hashlib.sha256(old_renderer.read_bytes()).hexdigest(),
       displayed_samples={'knee_angle_and_MN':int(keep.sum()),'knee_current':int(((DATA['knee']['t']>=29.85)&(DATA['knee']['t']<=30.25)).sum()),'cup':{cond:{key:len(DATA['cup'][k][key]['t'])for key in ['th','mCup','A','At']}for k,cond in enumerate(['on','off'])}},
       input_trace_npz=str(TRACE),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in source},
       knee_windows_s={'current':[29.85,30.25],'angle_and_MN':[29,40]},cup_window_s=[0,10],
       transformations=['native internal knee angle negated and converted radians to degrees','cup positive sag negated and converted radians to degrees for elevation','knee time offset minus30s for pulse clock','BPA glyph coordinates use exact sns_icon_geom.m expressions'],
       smoothing='none',phase_shift='none; only stated pulse-relative clock origin',decimation='none',interpolation='none beyond the source saved1ms resampling',path_simplify=False,
       output_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in OUT.glob('sns_*')if p.is_file()})
    (HERE/'CAD_SNS_FIGURE_PROVENANCE.json').write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf-8')
