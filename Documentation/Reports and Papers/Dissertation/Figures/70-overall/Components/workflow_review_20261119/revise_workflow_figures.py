"""Publish the dissertation workflow figures with arrows specified in print points.

Run from any directory with PyMuPDF and the Windows Arial fonts installed.
Input PDFs are the archived pre-revision versions in each figure directory.
"""
from pathlib import Path
import math
import os
import pymupdf as pdf

FIGURES = Path(__file__).resolve().parents[3]
FONTS = Path(os.environ.get("WINDIR", r"C:\Windows")) / "Fonts"
ARIAL = FONTS / "arial.ttf"
BOLD = FONTS / "arialbd.ttf"
INK = (0.13, 0.20, 0.25)
MUTED = (0.25, 0.32, 0.36)
BLUE = (0.0, 0.38, 0.59)
PALE = (0.95, 0.97, 0.98)
WHITE = (1, 1, 1)
GREEN = (23/255, 77/255, 36/255)
APP_BLUE = (66/255, 88/255, 105/255)
HEAD_LENGTH = 10.0
HEAD_WIDTH = 8.0
STROKE = 1.6

fonts = {False: pdf.Font(fontfile=str(ARIAL)), True: pdf.Font(fontfile=str(BOLD))}
files = {False: ARIAL, True: BOLD}


def text(pg, x, y, label, size=11, bold=False, align="center", color=INK):
    font = fonts[bold]
    name = "ArialBold" if bold else "ArialText"
    for i, line in enumerate(label.split("\n")):
        width = font.text_length(line, fontsize=size)
        start = x - (width/2 if align == "center" else width if align == "right" else 0)
        pg.insert_text((start, y+i*size*1.15), line, fontsize=size, fontname=name,
                       fontfile=str(files[bold]), color=color)


def box(pg, xy, title, body=None, future=False, size=11):
    x, y, w, h = xy
    pg.draw_rect(pdf.Rect(x,y,x+w,y+h), color=INK, fill=PALE, width=1.0,
                 dashes="[4 3] 0" if future else None)
    lines = [] if body is None else body.split("\n")
    block = (len(title.split("\n")) + len(lines))*size*1.15 + (3 if lines else 0)
    baseline = y+(h-block)/2+size
    text(pg,x+w/2,baseline,title,size,bold=bool(lines))
    if lines:
        text(pg,x+w/2,baseline+len(title.split("\n"))*size*1.15+3,body,size)


def line(pg, points, color=INK, dashed=False):
    pg.draw_polyline(points, color=color, width=STROKE, lineJoin=1,
                     dashes="[4 3] 0" if dashed else None)


def arrow(pg, points, color=INK, dashed=False, length=HEAD_LENGTH, width=HEAD_WIDTH):
    """Use a filled tip and stop the shaft within it. Keep all turns outside nodes."""
    a, b = points[-2:]
    dx, dy = b[0]-a[0], b[1]-a[1]
    last = math.hypot(dx,dy)
    if last < length:
        raise ValueError(f"Short final arrow segment {points}")
    ux, uy = dx/last, dy/last
    base = (b[0]-length*ux,b[1]-length*uy)
    inside = (b[0]-(length-.8)*ux,b[1]-(length-.8)*uy)
    line(pg,[*points[:-1],inside],color,dashed)
    pg.draw_polyline([b,(base[0]-.5*width*uy,base[1]+.5*width*ux),
                        (base[0]+.5*width*uy,base[1]-.5*width*ux),b],
                     color=None,fill=color,closePath=True)


def save(doc, folder, name, alt=None):
    pg=doc[0]
    assert len(doc)==1 and pg.rect.width>0
    doc.subset_fonts()
    data=doc.tobytes(garbage=4,deflate=True,no_new_id=True)
    path=FIGURES/folder/(name+".pdf")
    if not path.exists() or path.read_bytes()!=data:
        tmp=path.with_suffix(".new")
        tmp.write_bytes(data)
        tmp.replace(path)
    if alt:
        path.with_suffix(".alt.txt").write_text(alt+"\n",encoding="utf-8")
    doc.close()


def design():
    doc=pdf.open(); pg=doc.new_page(width=468,height=322)
    box(pg,(10,11,135,54),"Human benchmark","Torque and range\nof motion (Ch. 2)")
    box(pg,(211,11,150,54),"Select actuators","Diameter, resting length,\nand actuator count")
    box(pg,(10,105,135,54),"BPA force model","Actuator forces\n(Ch. 3)")
    box(pg,(211,105,150,54),"Design routes","Corrected torque model\n(Ch. 4)")
    box(pg,(10,210,145,51),"Refine design","Redesign routes or\nreselect actuators")
    box(pg,(349,206,109,59),"Closed-loop\nverification","(Ch. 6)",future=True)
    box(pg,(343,290,115,27),"Commit to hardware",future=True,size=10.9)
    pg.draw_polyline([(210,235),(265,205),(320,235),(265,265),(210,235)],
                     color=INK,fill=WHITE,width=1.0,closePath=True)
    text(pg,265,232,"Meets\ntargets?",11.2)
    arrow(pg,[(145,38),(211,38)]);text(pg,176,29,"targets",10.7)
    arrow(pg,[(286,65),(286,105)])
    arrow(pg,[(145,131),(211,131)])
    arrow(pg,[(265,159),(265,205)])
    text(pg,276,183,"Predicted torque and motion",10.7,align="left")
    arrow(pg,[(210,235),(155,235)]);text(pg,181,225,"no",10.7)
    arrow(pg,[(320,235),(349,235)]);text(pg,335,225,"yes",10.7)
    arrow(pg,[(403,265),(403,290)])
    # The return path rejoins the actuator-selection input. The interrupted
    # vertical segment crosses the separate force-model signal without joining it.
    line(pg,[(82,210),(82,181),(181,181),(181,137)])
    line(pg,[(181,125),(181,38)])
    pg.draw_circle((181,38),2.0,color=None,fill=INK)
    save(doc,"70-overall","design_workflow",
         "Human torque and motion targets determine actuator selection. The BPA force model and actuator selection inform routing. A failed torque or motion check returns to redesign and actuator selection. Passing leads to proposed closed-loop verification and a hardware commitment. The dashed boxes identify the last two steps.")


def hierarchy():
    doc=pdf.open(); pg=doc.new_page(width=468,height=336)
    text(pg,9,16,"Proposed supraspinal additions",12.1,True,"left")
    box(pg,(9,30,134,55),"Basal ganglia","Action and mode\nselection",future=True)
    box(pg,(161,30,137,55),"Cerebellum","Adaptive gait gating\nand phase tuning",future=True)
    box(pg,(317,30,124,55),"Vestibular / ocular","Orientation sensing\nand visual correction",future=True)
    text(pg,9,137,"Spinal network (Ch. 6)",12.1,True,"left")
    box(pg,(9,148,183,54),"Two-level CPG","Rhythm and pattern formation")
    box(pg,(243,148,185,54),"Sensorimotor layer","Motor pools; Ia, Ib, II;\nheel and toe contact")
    text(pg,9,257,"Simulated mechanical model",12.1,True,"left")
    box(pg,(9,271,419,31),"BPA actuators and robot model")
    arrow(pg,[(76,85),(76,99),(144,99),(144,148)],dashed=True)
    arrow(pg,[(230,85),(230,112),(175,112),(175,148)],dashed=True)
    arrow(pg,[(379,85),(379,148)],dashed=True)
    text(pg,273,117,"Descending\ncommands",10.8)
    arrow(pg,[(192,166),(243,166)])
    arrow(pg,[(243,190),(218,190),(218,234),(99,234),(99,202)])
    text(pg,154,225,"Sensory feedback",10.8)
    arrow(pg,[(279,202),(279,271)])
    text(pg,271,244,"Activation",10.8,align="right")
    arrow(pg,[(346,271),(346,202)])
    text(pg,354,230,"Proprioceptive\nand contact input",10.8,align="left")
    # A future sensor measures body orientation. Enter it from the outside.
    arrow(pg,[(428,286),(456,286),(456,58),(441,58)],MUTED,dashed=True)
    text(pg,10,324,"Balance-platform hardware is a precursor to orientation feedback.",10.8,align="left",color=MUTED)
    save(doc,"70-overall","sensorium_hierarchy",
         "Dashed future inputs lead from basal-ganglia and cerebellar modules to the two-level CPG, and from vestibular and ocular sensing to the sensorimotor layer. A dashed outer line returns body orientation to the proposed sensors. Solid arrows connect the CPG and sensorimotor layers in both directions and carry activation to the simulated BPA model and proprioceptive and contact input back to the sensorimotor layer. Earlier balance-platform hardware provided an orientation-feedback precedent.")


def pipeline():
    doc=pdf.open(); pg=doc.new_page(width=468,height=275)
    for x,title in [(76,"Model conversion"),(236,"Closed-loop simulation"),(395,"Future hardware")]:
        text(pg,x,17,title,12,True)
    pg.draw_rect(pdf.Rect(332,29,461,270),color=MUTED,width=.9,dashes="[4 3] 0")
    box(pg,(10,37,125,53),"OpenSim model","Gait2392 routing")
    box(pg,(10,119,125,35),"MyoConverter")
    box(pg,(10,192,125,56),"MuJoCo model","MJCF; activation\nstate actuators")
    box(pg,(181,37,111,59),"MuJoCo physics","BPA force model\n(Ch. 3)")
    box(pg,(181,207,111,46),"SNS-Toolbox","Neural controller")
    box(pg,(343,37,97,52),"Jetson Orin","SNS evaluation",future=True)
    box(pg,(343,130,97,54),"Joint controllers","Valves, pressure,\nand contact",future=True)
    box(pg,(343,228,97,31),"Treadmill biped",future=True)
    arrow(pg,[(73,90),(73,119)])
    arrow(pg,[(73,154),(73,192)])
    arrow(pg,[(135,220),(150,220),(150,67),(181,67)])
    arrow(pg,[(181,232),(165,232),(165,112),(208,112),(208,96)],BLUE)
    text(pg,185,136,"Per-muscle\nactivation",10.8,align="left",color=BLUE)
    arrow(pg,[(271,96),(271,112),(309,112),(309,186),(268,186),(268,207)])
    text(pg,185,166,"Ia, Ib, II, and\ncontact feedback",10.8,align="left")
    arrow(pg,[(292,232),(322,232),(322,53),(343,53)],dashed=True)
    arrow(pg,[(392,89),(392,130)])
    text(pg,379,105,"SNS\noutput",10.5,align="right")
    arrow(pg,[(392,184),(392,228)])
    text(pg,379,202,"Valve\ndrive",10.5,align="right")
    arrow(pg,[(440,244),(454,244),(454,160),(440,160)],MUTED)
    arrow(pg,[(440,144),(454,144),(454,73),(440,73)],MUTED,dashed=True)
    save(doc,"70-overall","pipeline",
         "OpenSim, MyoConverter, and the prepared MuJoCo model lead into a closed loop. SNS-Toolbox sends per-muscle activation up to MuJoCo physics and receives Ia, Ib, II, and contact feedback. The future hardware path connects the controller to a Jetson Orin, the joint controllers, and the treadmill biped. Distinct return paths run from the biped to joint controllers and from the controllers to the Orin.")


def replace_gap(pg,rect,points,color):
    pg.draw_rect(pdf.Rect(rect),color=None,fill=WHITE)
    arrow(pg,points,color)


def overlay(folder,name):
    src=FIGURES/folder/"Deprecated"/"20261119_toolchain"/"active_before"/(name+".pdf")
    doc=pdf.open(src); pg=doc[0]
    if name=="simulation_workflow_authoritative":
        for y in (42.2,120.1):
            for left,right in ((109.5,128.2),(318.6,337.7)):
                replace_gap(pg,(left-.3,y-5.4,right+.4,y+5.4),[(left,y),(right,y)],GREEN)
    elif name=="controller_retest_and_optimization":
        # Original circuit export. Keep each synapse block and neuron symbol,
        # and make the four signal directions as readable as the workflow.
        neuron_ink = (36/255,54/255,66/255)
        for y in (70.1,136.1):
            for tip in (88.5,152.2):
                pg.draw_rect(pdf.Rect(tip-11,y-4.2,tip+.15,y+4.2),color=None,fill=WHITE)
                arrow(pg,[(tip-14,y),(tip,y)],neuron_ink)
        for start,end in ((77.6,90.8),(145.1,159.9)):
            replace_gap(pg,(344,start-.10,356,end+.15),[(350,start),(350,end)],GREEN)
        for left,right in ((143.5,160.1),(305.9,322.5)):
            replace_gap(pg,(left-.1,267,right+.1,278),[(left,272.6),(right,272.6)],GREEN)
    elif name=="research_toolchain_workflow":
        # Preserve the six original boxes and the proposed interface band.
        for l,r,y in ((134.5,162.9,85.3),(295.1,323.4,85.3),
                      (323.4,295.1,177.5),(162.9,134.5,177.5)):
            a,b=sorted((l,r))
            replace_gap(pg,(a+.75,y-6,b-.75,y+6),
                        [(l+(1 if r>l else -1),y),(r-(1.1 if r>l else -1.1),y)],APP_BLUE)
        replace_gap(pg,(382,118.6,397,144),[(389.5,119.0),(389.5,143.7)],APP_BLUE)
    else:
        raise ValueError(name)
    save(doc,folder,name)


if __name__=="__main__":
    design()
    hierarchy()
    pipeline()
    overlay("60-sns","simulation_workflow_authoritative")
    overlay("60-sns","controller_retest_and_optimization")
    overlay("95-AppendixD","research_toolchain_workflow")
