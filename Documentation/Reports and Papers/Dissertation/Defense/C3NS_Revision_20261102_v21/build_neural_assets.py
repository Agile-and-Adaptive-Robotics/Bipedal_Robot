"""Rearrange saved vector plots for slide readability; no model execution."""
from pathlib import Path
from copy import deepcopy
import hashlib, json, re
import xml.etree.ElementTree as ET
import fitz

BASE = Path(__file__).resolve().parent / 'assets' / 'neural'
OUT = BASE / 'slide_panels'
OUT.mkdir(exist_ok=True)
S = 'http://www.w3.org/2000/svg'
ET.register_namespace('', S)
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
q = lambda name: '{'+S+'}'+name
manifest=[]

def setup(source, width, height):
    root=ET.parse(source).getroot()
    dst=ET.Element(q('svg'), {'width':str(width)+'pt','height':str(height)+'pt','viewBox':f'0 0 {width} {height}','version':'1.1'})
    defs=ET.SubElement(dst,q('defs'))
    seen=set()
    for d in root.iter(q('defs')):
        for elem in d:
            key=elem.get('id') or ET.tostring(elem)
            if key not in seen: defs.append(deepcopy(elem)); seen.add(key)
    ET.SubElement(dst,q('rect'), {'x':'0','y':'0','width':str(width),'height':str(height),'fill':'white'})
    return root,dst

def label(dst,text,x,y,size=14,bold=False,anchor='start'):
    t=ET.SubElement(dst,q('text'),{'x':str(x),'y':str(y),'font-family':'Arial','font-size':str(size),'font-weight':'bold' if bold else 'normal','text-anchor':anchor,'fill':'#202A25'})
    t.text=text

def bounds(axis):
    path=next(c for c in axis if c.get('id','').startswith('patch_')).find(q('path'))
    nums=list(map(float,re.findall(r'-?\d+(?:\.\d+)?',path.get('d'))))
    xs,ys=nums[0::2],nums[1::2]
    return min(xs),min(ys),max(xs),max(ys)

def find(root, id):
    return next(e for e in root.iter() if e.get('id')==id)

def plot(root,dst,axis_id,x,y,font=14,xticks=True,title=False):
    src=find(root,axis_id); a=deepcopy(src)
    xmin,ymin,xmax,ymax=bounds(src)
    for parent in list(a.iter()):
        for child in list(parent):
            if child.tag==q('defs'): parent.remove(child)
    # Replace only x-axis labels/grid; actual plotted paths and y scales stay exact.
    axes=[c for c in a if c.get('id','').startswith('matplotlib.axis_')]
    if axes: a.remove(axes[0])
    for t in a.iter(q('text')):
        t.set('style',re.sub(r'font:\s*[\d.]+px',f'font: {font}px',t.get('style','')))
    if not title:
        for child in list(a):
            if child.get('id','').startswith('text_'): a.remove(child)
    wrap=ET.SubElement(dst,q('g'), {'transform':f'translate({x-xmin} {y-ymin})'})
    wrap.append(a)
    if xticks:
        # Tick positions are copied from the source's 8/10/12-second data scale.
        for value in (8,10,12):
            tx=x+(value-7.663)/(13-7.663)*(xmax-xmin)
            # Neural source limits are exactly those of its first/right axis.
            alltext=[t for t in root.iter(q('text')) if t.text==str(value) and xmin <= float(t.get('x','-999')) <= xmax]
            if alltext:
                source_x=float(alltext[0].get('x'))
                if xmin <= source_x <= xmax: tx=x+source_x-xmin
            ET.SubElement(dst,q('line'), {'x1':str(tx),'x2':str(tx),'y1':str(y+ymax-ymin),'y2':str(y+ymax-ymin+3),'stroke':'black','stroke-width':'.8'})
            label(dst,str(value),tx,y+ymax-ymin+18,font,anchor='middle')
        label(dst,'Recorded time (s)',x+(xmax-xmin)/2,y+ymax-ymin+38,font,anchor='middle')
    return xmax-xmin,ymax-ymin

def save(dst,name,source,scope):
    file=OUT/(name+'.svg'); ET.ElementTree(dst).write(file,encoding='utf-8',xml_declaration=True)
    # SVG is the editable vector deliverable. Raster previews use Sharp so dashed
    # line styles survive; MuPDF's SVG reader drops dash patterns here.
    manifest.append({'file':str(file),'source':str(source),'source_sha256':hashlib.sha256(Path(source).read_bytes()).hexdigest(),'sha256':hashlib.sha256(file.read_bytes()).hexdigest(),'scope':scope,'data_changes':'none; retained vector data paths and y scales, translated axes, enlarged labels, restored shared time ticks'})

source=BASE/'ground9_rc_s4/PRESENTATION_NEURAL.svg'
root,dst=setup(source,750,245)
panels=[('axes_2','Rhythm generation',[('Extensor','#0072b2',False),('Flexor','#d55e00',False)]),('axes_5','RF/HAM pattern formation',[('RF PF','#0072b2',False),('HAM PF','#d55e00',False),('RF relay','#009e73',True),('HAM relay','#cc79a7',True)]),('axes_11','Knee motor neurons',[('Vastus intermedius','#0072b2',False),('Biceps femoris SH','#d55e00',False)])]
for i,(axis_id,title,legend) in enumerate(panels):
    px=i*250
    label(dst,title,px+143,23,16,True,'middle')
    for j,(name,color,dash) in enumerate(legend):
        lx=px+45+(j%2)*104 if len(legend)>2 else px+45
        ly=43+(j//2)*17 if len(legend)>2 else 43+j*17
        ln=ET.SubElement(dst,q('line'),{'x1':str(lx),'x2':str(lx+16),'y1':str(ly-4),'y2':str(ly-4),'stroke':color,'stroke-width':'1.4'})
        if dash:ln.set('stroke-dasharray','5,3')
        label(dst,name,lx+21,ly,12.5)
    plot(root,dst,axis_id,px+49.14,90,font=14.5)
save(dst,'NEURAL_RIGHT_HIERARCHY',source,'Right leg only, selected ground9 RC S4; RG, RF/HAM PF, and knee MN from the same recorded three cycles')

source=BASE/'contact_feedback/STAGE5_CONTACT_STRENGTH.svg'
root,dst=setup(source,468,274)
label(dst,'Quarter strength',146,19,14,True,'middle')
label(dst,'Original strength',365,19,14,True,'middle')
for axis_id,x,y in [('axes_3',61,58),('axes_4',278,58),('axes_13',61,164),('axes_14',278,164)]:
    w,h=plot(root,dst,axis_id,x,y,font=11.5,xticks=axis_id in ('axes_13','axes_14'),title=False)
    label(dst,'Ankle PF activity' if axis_id in ('axes_3','axes_4') else 'Vertical foot force',x, y-8,12)
save(dst,'CONTACT_PF_AND_FOOT_LOADING',source,'Matched quarter-strength versus original-strength S5; same controller/body/support settings apart from contact-PF scaling; original PF/DF and right/left color legends retained')

# Retain vector geometry from an existing PDF when isolating its final force panel.
source=BASE/'ground9_rc_s4/ALIGNED_MECHANICS.pdf'
doc=fitz.open(source); page=doc[0]
blocks=page.get_text('blocks')
title=next(b for b in blocks if 'Foot loads and assistance' in b[4])
clip=fitz.Rect(page.rect.width/2, title[1]-4, page.rect.width, page.rect.height)
out=fitz.open(); p=out.new_page(width=clip.width,height=clip.height)
p.show_pdf_page(p.rect,doc,0,clip=clip)
out.save(OUT/'FOOT_LOADS_AND_ASSISTANCE.pdf')
p.get_pixmap(matrix=fitz.Matrix(3,3),alpha=False).save(OUT/'FOOT_LOADS_AND_ASSISTANCE.png')
manifest.append({'file':str(OUT/'FOOT_LOADS_AND_ASSISTANCE.pdf'),'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'clip':list(clip),'scope':'Exact bottom-right source PDF panel including legend, units, and time axis; no chart/data editing'})
(OUT/'DERIVATIVE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')

# Source-preserving close-ups of the newly delivered motor-loop diagrams.
# They are explicitly labeled local crops of the historical standing export.
source=BASE/'eb475ws4_motor_loops/biarticular_full_loops.svg'
if source.exists():
    for name,viewbox in [('RF_MOTOR_FORCE_CROP',(70,1110,880,455)),('RF_LENGTH_SENSORY_CROP',(70,405,880,445))]:
        root=ET.parse(source).getroot()
        background=root.find(q('rect'))
        if background is not None:
            background.set('width','2400');background.set('height','1900')
        root.set('viewBox',' '.join(map(str,viewbox)))
        root.set('width',str(viewbox[2]));root.set('height',str(viewbox[3]))
        target=OUT/(name+'.svg')
        ET.ElementTree(root).write(target,encoding='utf-8',xml_declaration=True)
        manifest.append({'file':str(target),'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'viewBox':viewbox,'scope':'Local RF crop of saved w2lvar standing export, not current ground9 controller; paths extending beyond the viewport are outside the displayed subnetwork','data_changes':'none; original vector shapes and text unchanged'})
(OUT/'DERIVATIVE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps([{'asset':Path(v['file']).name,'scope':v['scope']} for v in manifest],indent=2))
