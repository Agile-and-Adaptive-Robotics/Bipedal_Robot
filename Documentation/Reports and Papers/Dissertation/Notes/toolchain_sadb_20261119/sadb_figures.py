"""Typeset the Chapter 5 circuit rule and enlarge selected portions of actual app captures."""
from pathlib import Path
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

here = Path(__file__).resolve().parent
folder = here.parents[1] / "Figures/50-database"
scratch = here.parents[4] / "tmp/toolchain_sadb_20261119"
captures = folder / "components/sadb_app_20261119"
plt.rcParams.update({"font.family":"Arial", "pdf.fonttype":42, "pdf.compression":9})
ink,blue,muted = "#182733","#204c72","#42545f"
width = 468

def page(height):
    f=plt.figure(figsize=(width/72,height/72),facecolor="white")
    ax=f.add_axes([0,0,1,1],xlim=(0,width),ylim=(height,0),aspect="equal")
    ax.set_axis_off()
    return f,ax

def label(ax,x,y,value,size=12,**kwargs):
    opts=dict(ha="left",va="center",fontsize=size,color=ink,linespacing=1.18)
    opts.update(kwargs)
    return ax.text(x,y,value,**opts)

def crop(ax,img,rect,x,y,scale):
    piece=img.crop(rect)
    w,h=piece.size
    ax.imshow(piece,extent=[x,x+w*scale,y+h*scale,y],interpolation="none",
              origin="upper",aspect="equal",zorder=1)

def save(f,stem,alt):
    f.savefig(folder/(stem+".pdf"),metadata={"Creator":"Chapter 5 figure generator","CreationDate":None})
    (folder/(stem+".alt.txt")).write_text(alt+"\n",encoding="utf-8")
    scratch.mkdir(parents=True,exist_ok=True)
    f.savefig(scratch/(stem+".png"),dpi=130)
    plt.close(f)

def record():
    grid=Image.open(captures/"pathway_table.png")
    drawer=Image.open(captures/"mccrea_record_identification.png")
    f,ax=page(337)
    label(ax,6,12,"A   Choose a pathway",size=12.5,fontweight="bold")
    crop(ax,grid,(987,21,1530,85),6,26,.465)
    label(ax,294,12,"B   Select the paper",size=12.5,fontweight="bold")
    crop(ax,grid,(21,255,218,438),294,26,.45)
    ax.plot([6,460],[116,116],color="#b4c0c7",lw=.8)
    label(ax,6,131,"C   Read the study and check its classifications",size=12.5,fontweight="bold")
    crop(ax,drawer,(21,7,753,194),6,146,.47)
    crop(ax,drawer,(23,301,443,384),6,254,.47)
    crop(ax,drawer,(23,725,558,907),210,243,.47)
    save(f,"sadb_table_record_guide","Actual SADb Table filter for Ib disynaptic excitation. Six entries match; selecting McCrea opens the paper’s identification. Enlarged portions of the drawer show the title, authors, DOI, the Ia and Ib afferent tags, and both Ib disynaptic inhibition and excitation pathway tags.")

def pivot():
    screen=Image.open(captures/"ib_cat_pivot_full.png")
    f,ax=page(319)
    label(ax,6,12,"A   Choose the row and column variables",size=12.5,fontweight="bold")
    crop(ax,screen,(25,10,901,79),6,24,.49)
    label(ax,6,73,"B   Select the cell to open the matching papers",size=12.5,fontweight="bold")
    cell=Image.open(captures/"ib_cat_pivot_cell.png")
    crop(ax,cell,(0,0,668,326),6,88,.477)
    label(ax,6,259,"C   Read the records or return to the pivot",size=12.5,fontweight="bold")
    drill=Image.open(captures/"ib_cat_drillthrough.png")
    crop(ax,drill,(21,158,459,229),6,275,.46)
    crop(ax,drill,(21,268,472,316),231,281,.46)
    save(f,"sadb_pivot_guide","Actual Pivot controls with Afferent types as rows and Animal as columns. The highlighted Ib by Cat cell contains 21 Papers. Selecting it opens the contributing records with the Back to pivot control.")

def rule():
    f,ax=page(178)
    def box(x,y,w,h,words):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0,rounding_size=5",
                     edgecolor=muted,facecolor="#f0f4f5",linewidth=1.15))
        label(ax,x+w/2,y+h/2,words,12,ha="center")
    for ypos,title,action,positive in ((6,"At rest","Inhibition",False),(93,"During the extensor phase","Excitation",True)):
        label(ax,6,ypos+6,title,13,fontweight="bold")
        y=ypos+25
        box(7,y,125,43,"Plantaris nerve\ngroup I input")
        box(297,y,164,43,"Medial gastrocnemius\nmotor neurons")
        path_y=y+26
        label(ax,214,y+8,action,11.2,ha="center",color=blue)
        ax.plot([132,292 if positive else 288],[path_y,path_y],color=blue,lw=2,solid_capstyle="butt")
        if positive:
            ax.add_patch(Polygon([[296,path_y],[284,path_y-5],[284,path_y+5]],color=blue))
        else:
            ax.plot([288,288],[path_y-7,path_y+7],color=blue,lw=2.5,solid_capstyle="butt")
    ax.plot([7,460],[82,82],color="#b4c0c7",lw=.8)
    save(f,"sadb_group_i_rule","The same plantaris group I stimulus inhibited medial gastrocnemius motor neurons at rest, shown by a barred line, and excited them during the extensor phase of midbrain-evoked fictive locomotion, shown by an arrowhead. These signs represent the recorded net motor-neuron response.")

if __name__=="__main__":
    record()
    pivot()
    rule()
