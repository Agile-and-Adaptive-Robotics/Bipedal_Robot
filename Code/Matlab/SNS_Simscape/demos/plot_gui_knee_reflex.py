"""Plot the native paired knee trials and draw the corresponding actuators."""
from argparse import ArgumentParser
from pathlib import Path
import json
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

BLUE = "#4477AA"
RED = "#AA444F"
GRAY = "#444444"


def trace(path):
    d = loadmat(path, simplify_cells=True)
    return d["time"] - d["onset"], dict(zip(d["signals"],d["samples"].T))


def connection(ax, pts, excitatory, c=GRAY):
    xy=np.asarray(pts,dtype=float)
    ax.plot(*xy.T, ls="-" if excitatory else "--", lw=1.16, color=c,zorder=1)
    delta=xy[-1]-xy[-2];delta/=np.linalg.norm(delta)
    normal=np.array([-delta[1],delta[0]]);q=xy[-1]
    if excitatory:
        ax.add_patch(Polygon([q+.15*delta,q-.12*delta-.105*normal,q-.12*delta+.105*normal],
                             facecolor="white",edgecolor=c,lw=1.0,zorder=4))
    else:
        ax.add_patch(Circle(q,.093,color=c,zorder=4))


def muscle(ax,center_x,center_y,geometry,scale):
    def xy(x,y):return np.column_stack((center_x+scale*np.asarray(x),center_y+scale*np.asarray(y)))
    ax.add_patch(Polygon(xy(geometry["bodyX"],geometry["bodyY"]),facecolor=geometry["fill"],edgecolor="black",lw=.75))
    for x0,y0,x1,y1 in geometry.get("fibers",geometry.get("braid",[])):
        ax.plot(*xy([x0,x1],[y0,y1]).T,color="black",lw=.56)
    for field in ("tendonR","tendonL","capR","capL"):
        if field in geometry:
            a=np.asarray(geometry[field])
            ax.add_patch(Polygon(xy(a[:,0],a[:,1]),facecolor=geometry.get("tendon",geometry.get("cap")),edgecolor="black",lw=.6))


def diagram(ax, biological, geometry):
    ax.set_xlim(-.1,14.1);ax.set_ylim(-.05,5.05);ax.set_aspect("equal");ax.axis("off")
    for i,(x,y,text) in enumerate([(2.5,3.75,"Ia"),(.75,2.85,"Ib"),(11.5,3.75,"Ia"),(13.25,2.85,"Ib")]):
        ax.add_patch(Circle((x,y),.31,ec="black",fc="#F0F5EE",lw=1.05,zorder=3))
        ax.text(x,y,text,ha="center",va="center",weight="bold",fontsize=10 if biological else 9.8,zorder=5)
    motors=[3.2,10.8];my=1.80
    for x in motors:
        ax.add_patch(Circle((x,my),.40,ec="black",fc="white",lw=1.16,zorder=3))
        ax.text(x,my,"MN",ha="center",va="center",weight="bold",fontsize=10 if biological else 9.4,zorder=4)
    connection(ax,[(2.50,3.44),(3.08,2.19)],True,BLUE)
    connection(ax,[(1.06,2.85),(2.73,1.85)],False,BLUE)
    connection(ax,[(11.50,3.44),(10.93,2.19)],True,RED)
    connection(ax,[(12.94,2.85),(11.27,1.85)],False,RED)
    connection(ax,[(2.81,3.75),(6.35,3.75),(6.35,1.84),(10.33,1.84)],False,BLUE)
    connection(ax,[(11.19,3.75),(7.58,3.75),(7.58,2.13),(3.62,2.13)],False,RED)
    for x in motors:
        ax.plot([x,x],[1.39,.94],color="black",lw=.95);muscle(ax,x,.70,geometry,.70)
    ax.text(2.05,4.54,"Extensor",ha="center",weight="bold",color=BLUE,fontsize=10.4)
    ax.text(11.97,4.54,"Flexor",ha="center",weight="bold",color=RED,fontsize=10.4)
    if biological:
        labels=("Vastus medialis","Biceps femoris short head")
    else:
        labels=("20 mm extensor BPA","20 mm flexor BPA")
    for x,label in zip(motors,labels):ax.text(x,.09,label,ha="center",va="center",fontsize=10 if biological else 9.3)
    ax.text(7,4.49,"Solid: excitation   Dashed: inhibition",ha="center",fontsize=10 if biological else 8.8)


def make(family,folder,output,icons):
    biological=family=="biological"
    current="ia_pulse" if biological else "nonspiking_pulse"
    control="control" if biological else "nonspiking_control"
    t,p=trace(folder/(current+".mat"));tc,c=trace(folder/(control+".mat"));assert np.array_equal(t,tc)
    delta={k:p[k]-c[k] for k in p}
    plt.rcParams.update({"font.family":"Arial","font.size":10.,"axes.linewidth":1.16,"pdf.fonttype":42,
                         "xtick.major.width":1.16,"ytick.major.width":1.16,"path.simplify":False})
    w=14*72/2.54;h=477. if biological else 529.;fig=plt.figure(figsize=(w/72,h/72))
    def txt(x,y,value,**kw):fig.text(x/w,y/h,value,**kw)
    def box(x,y,dx,dy):return fig.add_axes([x/w,y/h,dx/w,dy/h])
    txt(4,h-15,"(A)",weight="bold",size=10.6);txt(27,h-15,"Antagonist feedback circuit",weight="bold",size=10.4)
    if biological:
        diagram_box = (8,340,380,116)
        column_headings = (122.5,318,320)
        legend_anchor = (199/w,363.5/h)
        legend_text_size = 10
        # Give the two signed motor-response scales more vertical room.
        axes_rows = [(266,30),(201,43),(133.5,47),(84,29),(36,28)]
    else:
        diagram_box = (8,389,380,121)
        column_headings = (98,288,374)
        legend_anchor = (.54,.672)
        legend_text_size = 9
        axes_rows = [(280-60*k,32) for k in range(5)]
    diagram(box(*diagram_box),biological,icons["BIO" if biological else "BPA"])
    left,right,top=column_headings
    txt(left,top,"Pulse onset",weight="bold",ha="center");txt(right,top,"Response and recovery",weight="bold",ha="center")
    handles=[Line2D([],[],color=BLUE,lw=1.65,label="Extensor"),Line2D([],[],color=RED,ls="--",lw=1.55,label="Flexor")]
    fig.legend(handles=handles,loc="center",bbox_to_anchor=legend_anchor,ncol=2,frameon=True,fontsize=legend_text_size,
               title="Legend",title_fontproperties={"family":"Arial","weight":"bold","size":legend_text_size},
               handlelength=1.55 if biological else 1.8,columnspacing=1.0 if biological else 1.5,
               borderpad=.24,labelspacing=.18)
    rows=[("B","Injected Ia current (nA)",[(p["stim"],GRAY,"-")]),
          ("C","Change in motor-neuron voltage (mV)",[(delta["V_MN_ext"],BLUE,"-"),(delta["V_MN_flex"],RED,"--")])]
    if biological:
        rows.append(("D","Change in muscle activation (percentage points)",[(100*delta["A_ext"],BLUE,"-"),(100*delta["A_flex"],RED,"--")]))
    else:
        rows.append(("D","Change in BPA pressure (kPa)",[(delta["P_ext"],BLUE,"-"),(delta["P_flex"],RED,"--")]))
    rows.extend([("E","Change in "+("muscle" if biological else "actuator")+" force (N)",[(delta["F_ext"],BLUE,"-"),(delta["F_flex"],RED,"--")]),
                 ("F","Change in knee extension (deg)",[(-np.rad2deg(delta["th"]),GRAY,"-")])])
    spans=[(-.024,.250),(-.5,12 if biological else 4.5)]
    for k,(letter,heading,ss) in enumerate(rows):
        base,axis_height=axes_rows[k]
        header_y=base+axis_height+6
        txt(4,header_y,f"({letter})",weight="bold",size=10.1);txt(27,header_y,heading,weight="bold",size=10.0)
        wide=(t>=spans[1][0])&(t<=spans[1][1]);ys=np.concatenate([s[0][wide] for s in ss]);lo=min(0,float(min(ys)));hi=max(0,float(max(ys)))
        vspan=max(hi-lo,.01);lims=(lo-.10*vspan,hi+.13*vspan)
        for column,(xa,xb) in enumerate(spans):
            width = 140 if biological and column==1 else (151,143)[column]
            ax=box((47,248)[column],base,width,axis_height)
            ax.spines[["top","right"]].set_visible(False); ax.set_xlim(xa,xb);ax.set_ylim(lims)
            ax.axvline(0,color="#777777",ls=":",lw=.9,zorder=0)
            if k:ax.axhline(0,color=".80",lw=.7,zorder=0)
            ix=(t>=xa-.00001)&(t<=xb+.00001)
            for s,color,style in ss:
                if k==0:ax.step(t[ix],s[ix],where="post",color=color,lw=1.45)
                else:ax.plot(t[ix],s[ix],ls=style,color=color,lw=1.55 if style=="-" else 1.35)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=2,min_n_ticks=2,steps=[1,2,2.5,5,10]))
            if k==1 and biological:
                ax.set_ylim(min(-.205,lo-.02*vspan),max(.345,hi+.03*vspan))
                ax.set_yticks([-.15,0,.25])
            if k==2 and biological:
                ax.set_ylim(min(-1.11,lo-.02*vspan),max(3.08,hi+.03*vspan))
                ax.set_yticks([-1,0,2])
            ticks=[0,.05,.10,.15,.20,.25] if column==0 else ([0,3,6,9,12] if biological else [0,1,2,3,4])
            ax.set_xticks(ticks)
            ax.set_xticklabels([str(round(v*(1000 if column==0 else 1))) for v in ticks] if k==4 else [])
            ax.tick_params(labelsize=10 if biological else 9.5,pad=2,length=3.5)
    txt(122.5 if biological else 116,8,"Time from pulse onset (ms)",weight="bold",ha="center",size=10 if biological else 9.8)
    txt(318 if biological else 311,8,"Time from pulse onset (s)",weight="bold",ha="center",size=10 if biological else 9.8)
    stem="sns_knee_reflex_authoritative" if biological else "knee_bpa_reflex";output.mkdir(parents=True,exist_ok=True)
    metadata={"Title":"Pneumatic knee reflex"}
    if biological: metadata={"Title":"Simulated biological knee reflex","CreationDate":None}
    fig.savefig(output/(stem+".pdf"),metadata=metadata)
    fig.savefig(output/(stem+".png"),dpi=150);plt.close(fig)
    extension=-np.rad2deg(delta["th"]);peak=np.argmax(extension)
    desc=(f"(A) Each Ia sensory neuron excites its corresponding extensor or flexor motor neuron and inhibits the antagonist. Ib input inhibits its corresponding motor neuron. "
          f"The two actuator drawings show {('pennate human muscles' if biological else '20 mm braided pneumatic muscles')}. "
          "(B) A 2 nA current pulse is applied to extensor Ia for 30 ms. (C) Extensor motor-neuron voltage rises and flexor motor-neuron voltage falls. "
          f"(D) {( 'Muscle activation' if biological else 'Actuator pressure')} and (E) force changes are shown for each side. "
          f"(F) Knee extension peaks at {extension[peak]:.3f} degrees at {t[peak]*1000:.0f} ms. Changes are relative to the matched unstimulated model. "
          "The time-series panels share time zero at pulse onset in the detailed and longer views.")
    (output/(stem+".alt.txt")).write_text(desc+"\n",encoding="utf-8")


if __name__=="__main__":
    q=ArgumentParser();q.add_argument("family",choices=["biological","pneumatic"]);q.add_argument("results",type=Path)
    q.add_argument("output",type=Path);q.add_argument("icon_geometry",type=Path);a=q.parse_args()
    make(a.family,a.results,a.output,json.loads(a.icon_geometry.read_text()))
