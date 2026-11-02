"""Static connection drawings. Sources are parsed as data, never imported or run."""
import csv,html,json,math,re
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as E
import pymupdf
from audit_saved_graphs import ROOT,OUT,table,digest

DEST=OUT/'drawings/v3';DEST.mkdir(parents=True,exist_ok=True)
BASE=ROOT/'Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families'
def text(x,y,s,size=20,anchor='middle',color='#111'):
    return ''.join(f'<text x="{x}" y="{y+i*size*1.2}" font-family="Arial" font-size="{size}" text-anchor="{anchor}" fill="{color}">{html.escape(t)}</text>' for i,t in enumerate(s.split('\n')))
def dist(a,b):return math.hypot(a[0]-b[0],a[1]-b[1])
def toward(a,b,d):
    length=dist(a,b);return (a[0]+(b[0]-a[0])*d/length,a[1]+(b[1]-a[1])*d/length) if length else a
def pt(p):return f'{p[0]:.2f},{p[1]:.2f}'
def intersect(a,b,c,d):
    u=(b[0]-a[0],b[1]-a[1]);v=(d[0]-c[0],d[1]-c[1]);det=u[0]*v[1]-u[1]*v[0]
    if abs(det)<1e-8:return None
    t=((c[0]-a[0])*v[1]-(c[1]-a[1])*v[0])/det;s=((c[0]-a[0])*u[1]-(c[1]-a[1])*u[0])/det
    if 0<t<1 and 0<s<1:
        p=(a[0]+t*u[0],a[1]+t*u[1])
        if min(dist(p,z) for z in (a,b,c,d))>17:return p
    return None

class Canvas:
    def __init__(self,title,h=1390):
        self.title=title;self.h=h;self.nodes={};self.groups=[];self.wires=[];self.records=[];self.labels=[];self.bridge_count=0
    def node(self,name,label,x,y,w=230,h=70,color='#FFF0B3',shape='box'):
        self.nodes[name]={'name':name,'label':label,'x':x,'y':y,'w':w,'h':h,'color':color,'shape':shape}
    def group(self,label,x,y,w,h):self.groups.append((label,x,y,w,h))
    def wire(self,points,edges,sign=None,color='#333',dashed=False,label=None):
        points=[tuple(p) for i,p in enumerate(points) if i==0 or tuple(p)!=tuple(points[i-1])]
        assert len(points)>1
        if sign is None:sign=edges[0].get('sign_category','exc')
        self.wires.append({'points':points,'edges':edges,'sign':sign,'color':color,'dashed':dashed,'label':label})
        self.records.extend(edges)
    def output(self,slug,note=''):
        svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="1800" height="{self.h}" viewBox="0 0 1800 {self.h}"><rect width="100%" height="100%" fill="white"/>'+text(900,45,self.title,29)
        for label,x,y,w,h in self.groups:
            svg+=f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="16" fill="#F6FAF2" stroke="#333" stroke-width="1.6"/>'
            svg+=f'<rect x="{x+w*.22}" y="{y}" width="{w*.56}" height="36" rx="8" fill="#A7DF46" stroke="#333" stroke-width="1.4"/>'+text(x+w/2,y+25,label,22)
        previous=[];ends=[]
        for wire in self.wires:
            points=wire['points'];path='M '+pt(points[0]);radius=10
            for i,(a,b) in enumerate(zip(points,points[1:])):
                start=a if i==0 else toward(a,points[i-1],min(radius,dist(a,points[i-1])/3,dist(a,b)/3))
                finish=b if i==len(points)-2 else toward(b,a,min(radius,dist(a,b)/3,dist(b,points[i+2])/3))
                if i:
                    after=toward(a,b,min(radius,dist(a,b)/3,dist(a,points[i-1])/3));path+=' Q '+pt(a)+' '+pt(after);start=after
                jumps=[]
                for c,d in previous:
                    p=intersect(start,finish,c,d)
                    if p is not None and all(dist(p,q)>18 for q in jumps):jumps.append(p)
                jumps.sort(key=lambda p:dist(start,p))
                for p in jumps:
                    before=toward(p,start,7);after=toward(p,finish,7)
                    if min(dist(start,p),dist(finish,p))<18:continue
                    path+=' L '+pt(before)+' A 7 7 0 0 1 '+pt(after);self.bridge_count+=1
                path+=' L '+pt(finish)
            dash=' stroke-dasharray="9 6"' if wire['dashed'] else ''
            # White under-stroke makes every semicircular bridge a clean overpass.
            svg+=f'<path d="{path}" fill="none" stroke="white" stroke-width="6.5" stroke-linecap="round"/>'
            svg+=f'<path d="{path}" fill="none" stroke="{wire["color"]}" stroke-width="2" stroke-linejoin="round"{dash}/>'
            previous.extend(zip(points,points[1:]));ends.append(wire)
        for n in self.nodes.values():
            if n['shape']=='circle':svg+=f'<ellipse cx="{n["x"]+n["w"]/2}" cy="{n["y"]+n["h"]/2}" rx="{n["w"]/2}" ry="{n["h"]/2}" fill="{n["color"]}" stroke="#222" stroke-width="1.5"/>'
            else:svg+=f'<rect x="{n["x"]}" y="{n["y"]}" width="{n["w"]}" height="{n["h"]}" rx="12" fill="{n["color"]}" stroke="#222" stroke-width="1.5"/>'
            lines=len(n['label'].split('\n'));svg+=text(n['x']+n['w']/2,n['y']+n['h']/2-(lines-1)*12+7,n['label'],20)
        for wire in ends:
            end=wire['points'][-1];before=wire['points'][-2];tip=toward(end,before,15);length=dist(before,end);ux=(end[0]-before[0])/length;uy=(end[1]-before[1])/length
            if wire['sign']=='inh':
                center=toward(end,before,6);svg+=f'<circle cx="{center[0]}" cy="{center[1]}" r="6" fill="{wire["color"]}"/>'
            else:
                a=(end[0]-uy*7,end[1]+ux*7);b=(end[0]+uy*7,end[1]-ux*7)
                svg+=f'<polygon points="{pt(tip)} {pt(a)} {pt(b)}" fill="white" stroke="{wire["color"]}" stroke-width="1.6"/>'
            if wire['label']:
                p=wire['points'][-2];svg+=text(p[0],p[1]-12,wire['label'],17)
        for item in self.labels:svg+=text(*item)
        if note:svg+=text(900,self.h-45,note,18)
        svg+='</svg>'
        path=DEST/(slug+'.svg');path.write_text(svg,encoding='utf-8')
        source=pymupdf.open(stream=svg.encode(),filetype='svg');pdf=pymupdf.open('pdf',source.convert_to_pdf());pdf.save(DEST/(slug+'.pdf'));pdf[0].get_pixmap(matrix=pymupdf.Matrix(.85,.85)).save(DEST/(slug+'.png'));pdf.close();source.close()
        (DEST/(slug+'_connections.json')).write_text(json.dumps({'title':self.title,'source_connections':self.records,'semicircular_bridges':self.bridge_count,'circuits_executed':0},indent=2))
        return {'slug':slug,'connections':len(self.records),'bridges':self.bridge_count}

def native():
    p=ROOT/'tmp/functional_subnetwork_remote/W2L108_verified/editable/W2L_VirtualWalker_Newer_Standalone.asim';r=E.parse(p).getroot();ns={n.findtext('ID'):n for n in r.iter('Neuron')};types={t.findtext('ID'):t for t in r.iter('SynapseType')};es=[]
    for e in r.iter('Connexion'):
        a,b=e.findtext('SourceID'),e.findtext('TargetID');t=types[e.findtext('SynapseTypeID')];rev=float(t.findtext('Equil'))
        es.append({'id':e.findtext('ID'),'source':ns[a].findtext('Name'),'target':ns[b].findtext('Name'),'gmax_uS':float(t.findtext('SynAmp')),'reversal_mV':rev,'threshold_mV':float(t.findtext('ThreshV')),'saturation_mV':float(t.findtext('SaturateV')),'sign_category':'exc' if rev>float(ns[b].findtext('RestingPot')) else 'inh'})
    return p,es

def module(c,label,y,names,show_ins=False):
    c.group(label,630,y,540,220 if show_ins else 165)
    c.node(names[0],'ext HC',650,y+55,110,100,'#FFE52E');c.node(names[1],'flx HC',1040,y+55,110,100,'#FFE52E')
    if show_ins:
        c.node(names[2],'ext IN',790,y+165,90,42,'#D9EAC1');c.node(names[3],'flx IN',920,y+165,90,42,'#D9EAC1')

def draw_w2l():
    p,es=native();c=Canvas('W2L108 · grouped RG–PF feedback',1390)
    for label,y,stem in [('RG',160,'L RG'),('Hip PF',550,'L Hip PF'),('K&A PF',940,'L Knee PF')]:module(c,label,y,[stem+' ext',stem+' flx',stem+' ext IN',stem+' flx IN'],label=='RG')
    sensors={'L heel contact':(70,150),'L Hip ext Ib':(70,450),'L Knee ext Ib':(70,820),'L Hip flx Ia':(1480,440),'L Hip flx II':(1480,560),'L Knee flx Ia':(1480,830),'L Knee flx II':(1480,1000),'R c1':(1480,110),'R heel contact':(1480,260),'R V3':(850,65)}
    for name,(x,y) in sensors.items():
        label=name[2:] if name.startswith('L ') else 'Opposite '+name[2:]
        c.node(name,label,x,y,230 if name=='R V3' else 245,72)
    selected=[e for e in es if e['source'] in sensors and e['target'] in c.nodes]
    slots=defaultdict(list)
    for e in selected:slots[e['target']].append(e)
    per_source=defaultdict(int)
    for i,e in enumerate(selected):
        a=c.nodes[e['source']];b=c.nodes[e['target']];j=slots[e['target']].index(e);sy=a['y']+18+per_source[e['source']]*18;per_source[e['source']]+=1;ty=b['y']+18+j*16
        if e['source']=='R V3':
            points=[(a['x'],a['y']+34),(590,a['y']+34),(590,b['y']+20),(b['x'],b['y']+20)]
        elif e['source'] in ('R c1','R heel contact'):
            if a['y']+12<=ty<=a['y']+a['h']-12:points=[(a['x'],ty),(b['x']+b['w'],ty)]
            else:points=[(a['x'],a['y']+36),(1390,a['y']+36),(1390,ty),(b['x']+b['w'],ty)]
        elif a['x']<500:
            lane={'L heel contact':350,'L Hip ext Ib':460,'L Knee ext Ib':410}[e['source']]+per_source[e['source']]*25
            if a['y']+12<=ty<=a['y']+a['h']-12:points=[(a['x']+a['w'],ty),(b['x'],ty)]
            else:points=[(a['x']+a['w'],sy),(lane,sy),(lane,ty),(b['x'],ty)]
        else:
            lane={'L Hip flx Ia':1370,'L Hip flx II':1320,'L Knee flx Ia':1420,'L Knee flx II':1270}[e['source']]+per_source[e['source']]*22
            if a['y']+12<=ty<=a['y']+a['h']-12:points=[(a['x'],ty),(b['x']+b['w'],ty)]
            else:points=[(a['x'],sy),(lane,sy),(lane,ty),(b['x']+b['w'],ty)]
        c.wire(points,[e])
    # Anatomical group-to-group links bundle the two actual RG→PF projections.
    for stem,y,lane in [('L Hip PF',550,600),('L Knee PF',940,1210)]:
        feed=[e for e in es if e['source'] in ('L RG ext','L RG flx') and e['target'] in (stem+' ext',stem+' flx')]
        c.wire([(900 if y==550 else 950,380),(900 if y==550 else 950,440),(lane,440),(lane,y-45),(900,y-45),(900,y)],feed,color='#888',label='RG → PF ×2')
    c.labels.append((900,1295,'Excitation: triangle base on destination · Inhibition: filled circle · Crossings: semicircular bridges',18))
    return c.output('W2L108_feedback','Internal RG/PF wiring omitted. Opposite-side sources are placed beside their RG destinations.'),p,es

def load(family):
    p=BASE/(family+'_graph_ledger.json');d=json.loads(p.read_text());return p,d['chemical_edges'],d

def draw_variant():
    p,es,d=load('w2lvar');c=Canvas('w2lvar · RG, Hip PF, and K&A PF',1390)
    stems=[('RG',150,['RG_E_l','RG_F_l','InE_l','InF_l']),('Hip PF',550,['PF_HIP-E_l','PF_HIP-F_l','PF_IN_HIP-E_l','PF_IN_HIP-F_l']),('K&A PF',950,['PF_KNEE-E_l','PF_KNEE-F_l','PF_IN_KNEE-E_l','PF_IN_KNEE-F_l'])]
    for label,y,names in stems:module(c,label,y,names,True)
    c.node('Ib','Ib · 27 channels\nhip / knee / ankle / trunk',60,320,270,80);c.node('Ia','Ia · 19 channels\nhip / knee / ankle',1460,450,280,80);c.node('II','II · 19 channels\nhip / knee / ankle',1460,710,280,80)
    c.node('HEEL_l','Heel contact',60,105,270,80);c.node('HEEL_r','Opposite heel',1450,300,285,70)
    c.node('V0D_r','Opposite V0D',1450,95,285,70);c.node('V3E_r','Opposite V3',60,230,270,70)
    c.node('TOE_l','Toe contact',1440,1180,280,65);c.node('TOEDF_l','Toe relay\nTOEDF',1200,1180,180,65,'#D9EAC1')
    grouped=defaultdict(list)
    for e in es:
        if e['target'] not in c.nodes:continue
        src=e['source'];family=src.split('_',1)[0]
        if family in ('Ia','II','Ib') and src.endswith('_l'):grouped[family,e['target']].append(e)
        elif src in c.nodes and src in ('HEEL_l','HEEL_r','V0D_r','V3E_r','TOEDF_l','TOE_l'):grouped[src,e['target']].append(e)
    # Each source category gets distinct ports; every destination remains its real cell.
    order=defaultdict(int);target_order=defaultdict(int);totals=defaultdict(int);source_totals=defaultdict(int)
    for src,target in grouped:totals[target]+=1;source_totals[src]+=1
    for i,((src,target),edges) in enumerate(grouped.items()):
        a=c.nodes[src];b=c.nodes[target];k=order[src];order[src]+=1;t=target_order[target];target_order[target]+=1
        ty=b['y']+12+t*(b['h']-24)/max(1,totals[target]-1)
        if src=='TOE_l':points=[(1440,1212),(1380,1212)]
        elif src=='TOEDF_l':points=[(1200,1212),(1180,1192),(1180,ty),(b['x']+b['w'],ty)]
        elif src=='V0D_r':points=[(1450,130),(1370,130),(1370,ty),(1150,ty)]
        elif src=='HEEL_r':points=[(1450,345),(1010,345)]
        elif src=='V3E_r':
            sy=a['y']+18+k*22;lane=560+k*20;delta=abs(ty-sy);points=[(330,sy),(lane-delta,sy),(lane,ty),(b['x'],ty)]
        elif a['x']<500:
            sy=a['y']+12+k*(a['h']-24)/max(1,source_totals[src]-1);lane=(390 if src=='HEEL_l' else 480)+(source_totals[src]-1-k)*18
            points=[(a['x']+a['w'],sy),(lane,sy),(lane,ty),(b['x'],ty)]
        else:
            sy=a['y']+12+k*(a['h']-24)/max(1,source_totals[src]-1);lane=(1320 if src=='Ia' else 1220)+k*16
            points=[(a['x'],sy),(lane,sy),(lane,ty),(b['x']+b['w'],ty)]
        c.wire(points,edges,label=(src+' ×'+str(len(edges))) if len(edges)>1 else None)
    for target,source in [('PF_HIP-E_l','RG_E_l'),('PF_HIP-F_l','RG_F_l'),('PF_KNEE-E_l','RG_E_l'),('PF_KNEE-F_l','RG_F_l')]:
        edge=[e for e in es if e['source']==source and e['target']==target]
        a=c.nodes[source];b=c.nodes[target];x=580 if '-E_' in target else 1180
        c.wire([(a['x']+55,a['y']+a['h']),(a['x']+55,a['y']+a['h']+45),(x,a['y']+a['h']+45),(x,b['y']-35),(b['x']+55,b['y']-35),(b['x']+55,b['y'])],edge,color='#888')
    c.labels.append((900,1310,'×N = separate muscle-channel synapses; their individual destination parameters are retained in the key.',18))
    return c.output('w2lvar_feedback','K&A is the shared knee–ankle PF pair. This page uses one set of destination objects, without duplicate HCs.'),p,es

def draw_syn6():
    p,es,d=load('syn6');c=Canvas('syn6 · RG and six NMF pattern-formation channels',1390)
    module(c,'RG',140,['RG_E_l','RG_F_l','InE_l','InF_l'],True)
    c.group('Six NMF PF channels',580,540,640,690)
    positions={'PF_S2_l':(650,610),'PF_S3_l':(650,830),'PF_S6_l':(650,1050),'PF_S1_l':(1040,610),'PF_S4_l':(1040,830),'PF_S5_l':(1040,1050),'PF_IN_E_l':(810,590),'PF_IN_F_l':(810,1130)}
    for name,(x,y) in positions.items():c.node(name,('ext PF IN' if name=='PF_IN_E_l' else 'flx PF IN' if name=='PF_IN_F_l' else name.replace('_l','').replace('PF_','')),x,y,110,70,'#D9EAC1' if 'IN' in name else '#FFE52E')
    c.node('HEEL_l','Heel contact',65,330,270,75);c.node('V0D_r','Opposite V0D',1470,100,270,70);c.node('V3E_r','Opposite V3',65,100,270,70)
    c.node('TOE_l','Toe contact',1470,1050,270,70);c.node('TOEDF_l','Toe relay\nTOEDF',1270,1050,170,70,'#D9EAC1')
    selected=[e for e in es if e['source'] in ('RG_E_l','RG_F_l','HEEL_l','V0D_r','V3E_r','TOE_l','TOEDF_l','PF_IN_E_l','PF_IN_F_l') and e['target'] in c.nodes and (not e['source'].startswith('RG_') or e['target'].startswith('PF_'))]
    used=defaultdict(int)
    for e in selected:
        a=c.nodes[e['source']];b=c.nodes[e['target']];k=used[e['source']];used[e['source']]+=1;red=e['source']=='HEEL_l' and e['target']=='InF_l'
        if e['source']=='TOE_l':pts=[(1470,1085),(1440,1085)]
        elif e['source']=='TOEDF_l':pts=[(1270,1085),(1150,1085)]
        elif e['source']=='HEEL_l':
            if e['target']=='InE_l':pts=[(335,350),(500,350),(500,320),(790,320)]
            elif e['target']=='InF_l':pts=[(335,375),(520,375),(520,430),(965,430),(965,347)]
            else:pts=[(335,395),(470,395),(470,605),(810,605)]
        elif e['source']=='V0D_r':pts=[(1470,135),(1370,135),(1370,b['y']+30),(1150,b['y']+30)]
        elif e['source']=='V3E_r':
            sy=120+k*22;pts=[(335,sy),(560+k*28,sy),(560+k*28,b['y']+30),(b['x'],b['y']+30)]
        elif e['source'].startswith('RG_'):
            x=520+k*18 if e['source']=='RG_E_l' else 1260+k*18
            if e['source']=='RG_E_l' and e['target']=='PF_S5_l':pts=[(a['x']+55,a['y']+a['h']),(a['x']+55,390),(x,390),(x,1015),(900,1015),(900,1085),(1040,1085)]
            elif e['source']=='RG_F_l' and e['target']=='PF_S5_l':pts=[(a['x']+55,a['y']+a['h']),(a['x']+55,390),(x,390),(x,1000),(1095,1000),(1095,1050)]
            else:pts=[(a['x']+55,a['y']+a['h']),(a['x']+55,a['y']+a['h']+50+k*18),(x,a['y']+a['h']+50+k*18),(x,b['y']+35),(b['x'] if e['source']=='RG_E_l' else b['x']+b['w'],b['y']+35)]
        else:
            if e['source']=='PF_IN_F_l':pts=[(810,1150+k*15),(780,1150+k*15),(780,b['y']+35),(760,b['y']+35)]
            elif e['target']=='PF_S1_l':pts=[(920,645),(1040,645)]
            else:pts=[(920,615),(970,615),(970,b['y']+35),(1040,b['y']+35)]
        c.wire(pts,[e],color='#B13B3B' if red else '#888' if e['source'].startswith('RG_') else '#333',dashed=red)
    c.labels.append((900,1290,'Red dashed: the saved heel→flexor-IN branch disputed in Comment 30; the source circuit has not been changed.',18))
    return c.output('syn6_feedback','S1–S6 are PF channels, not hip/knee/ankle pairs. Saved contact-input gains are zero in this standing condition.'),p,es

def key_and_findings(native_edges,variant_edges,syn_edges):
    chosen=[e for e in variant_edges if e['source']=='II_tib_ant_l']
    local=[e for e in variant_edges if e['source'] in ('IIX_tib_ant_l','IIIN_tib_ant_l')]
    compare={ (e['source'],e['target']):e for e in variant_edges }
    fields=('gmax_uS','reversal_mV','threshold_mV','saturation_mV')
    rows=[]
    for e in chosen+local:
        mirror=compare.get((e['source'][:-2]+'_r',e['target'][:-2]+'_r'))
        rows.append({'source':e['source'],'destination':e['target'],**{k:e[k] for k in fields},'left_right_equal':mirror is not None and all(e[k]==mirror[k] for k in fields),'edge_id':e['id']})
    table(DEST/'II_tib_ant_destination_parameters.csv',rows)
    report={'comment14':'Opposite c1→RG flexor exists in saved W2L108 and is now included.','comments15_23_24':'Saved W2L108 heel targets extensor HCs, not extensor IN. Its ankle afferent edges target the local motor/reflex network. These are source discrepancies relative to Ben’s expected diagram; no ankle or heel→IN connection was silently invented.','comment28':'w2lvar PF_KNEE is explicitly the shared K&A PF pair per build_network_w2lvar.py; displayed as K&A.','comment29':'II_tib_ant names an afferent neuron. Each outgoing synapse is destination-specific. Four central edges currently share g=0.027894736842105264uS,Erev=8mV, while the local IIX/IIIN relay edges use1uS. Left/right equality is checked per table row, not assumed globally. Source IDs retain side suffixes; visual labels omit them.','comments30_31':'Saved syn6 contains HEEL_l→InE_l and HEEL_l→InF_l. The disputed latter edge is visibly marked red/dashed; no circuit or saved weights changed. No absent-route boxes are retained.','omitted_blank_page':'s3k absent-feedback page removed from this drawing pack.','II_tib_ant_rows':rows,'circuits_executed':0}
    (DEST/'SOURCE_FINDINGS_AND_COMMENT_RESPONSES.json').write_text(json.dumps(report,indent=2))
    # A concise textual key is kept separate from the diagrams.
    svg='<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="1400"><rect width="100%" height="100%" fill="white"/>'+text(750,45,'Source / destination / synapse: what the labels mean',28)
    statements=[(85,'II_tib_ant names one type-II afferent neuron, not one globally fixed synapse.'),(125,'A connection is defined by its source, destination, conductance, reversal, and presynaptic window.'),(165,'The same afferent can therefore have different parameters at different destinations.'),(205,'The checked left/right homologues below match; diagram labels omit the side suffix.')]
    for y,s in statements:svg+=text(50,y,s,22,'start')
    for s,x in [('Source → destination',45),('g (uS)',750),('Erev (mV)',910),('lo / hi (mV)',1090),('L = R',1350)]:svg+=text(x,285,s,20,'start')
    for i,row in enumerate(rows):
        y=330+i*48
        if i%2==0:svg+=f'<rect x="30" y="{y-28}" width="1440" height="46" rx="5" fill="#EFF3F4"/>'
        lab=lambda s:s.removesuffix('_l')
        svg+=text(45,y,lab(row['source'])+' → '+lab(row['destination']),19,'start')+text(750,y,f'{row["gmax_uS"]:.6g}',19,'start')+text(910,y,str(row['reversal_mV']),19,'start')+text(1090,y,f'{row["threshold_mV"]:g} / {row["saturation_mV"]:g}',19,'start')+text(1350,y,'yes' if row['left_right_equal'] else 'no',19,'start')
    svg+=text(50,1290,'A ×N line bundles N separately parameterized connections; it is not a single summed synapse.',22,'start')
    svg+=text(50,1330,'Ankle channels are included in the muscle bundles where present. Complete members and edge IDs accompany the drawings.',20,'start')+'</svg>'
    path=DEST/'synapse_key.svg';path.write_text(svg,encoding='utf-8');s=pymupdf.open(stream=svg.encode(),filetype='svg');d=pymupdf.open('pdf',s.convert_to_pdf());d.save(DEST/'synapse_key.pdf');d[0].get_pixmap(matrix=pymupdf.Matrix(.85,.85)).save(DEST/'synapse_key.png');d.close();s.close()

if __name__=='__main__':
    a,p,ne=draw_w2l();b,vp,ve=draw_variant();c,sp,se=draw_syn6();key_and_findings(ne,ve,se)
    pages=['W2L108_feedback','w2lvar_feedback','syn6_feedback','synapse_key'];d=pymupdf.open()
    for name in pages:
        with pymupdf.open(DEST/(name+'.pdf')) as source:d.insert_pdf(source)
    d.save(DEST/'RG_PF_Feedback_v3.pdf');d.close()
    (DEST/'DRAWING_MANIFEST.json').write_text(json.dumps({'drawings':[a,b,c],'sources':{str(x):digest(x) for x in (p,vp,sp)},'pages':pages,'circuits_executed':0},indent=2))
    print(json.dumps({'pages':4,'drawings':[a,b,c],'circuits_executed':0}))
