"""Grouped schematics and pathway views of saved data. No circuit execution."""
import copy,json,re
from pathlib import Path
import xml.etree.ElementTree as E
import pymupdf
import draw_feedback_v3 as V
from audit_saved_graphs import ROOT,OUT,digest

DEST=OUT/'drawings/v4';DEST.mkdir(parents=True,exist_ok=True)
V.DEST=DEST
ORIGINAL=V.Canvas
BASE_OUTPUT=ORIGINAL.output
def boundary(p,n):
    x,y,w,h=n['x'],n['y'],n['w'],n['h'];px,py=p
    return ((abs(px-x)<.1 or abs(px-x-w)<.1) and y-.1<=py<=y+h+.1) or ((abs(py-y)<.1 or abs(py-y-h)<.1) and x-.1<=px<=x+w+.1)
def checked_output(self,slug,note=''):
    objects=list(self.nodes.values())+[dict(x=x,y=y,w=w,h=h) for title,x,y,w,h in self.groups]
    dangling=[];obscured=[]
    for i,wire in enumerate(self.wires):
        points=wire['points']
        for endpoint in (points[0],points[-1]):
            if not any(boundary(endpoint,n) for n in objects):dangling.append({'wire':i,'point':endpoint})
        for n in self.nodes.values():
            for a,b in zip(points,points[1:]):
                for k in range(1,40):
                    p=(a[0]+(b[0]-a[0])*k/40,a[1]+(b[1]-a[1])*k/40)
                    if n['x']+3<p[0]<n['x']+n['w']-3 and n['y']+3<p[1]<n['y']+n['h']-3:
                        obscured.append({'wire':i,'node':n['name']});break
                else:continue
                break
    result={'dangling_endpoints':dangling,'wires_through_node_bodies':obscured,'wire_count':len(self.wires),'circuits_executed':0}
    (DEST/(slug+'_geometry_check.json')).write_text(json.dumps(result,indent=2))
    return BASE_OUTPUT(self,slug,note)
ORIGINAL.output=checked_output
CAPTURED=[]
class CaptureCanvas(ORIGINAL):
    def output(self,slug,note=''):
        record=super().output(slug,note);CAPTURED.append((slug,self,record));return record
V.Canvas=CaptureCanvas

def focus(original,slug,kind):
    c=copy.deepcopy(original)
    def keep(w):
        if w['color']=='#888':return True
        sources=[e['source'] for e in w['edges']]
        return any(s.startswith(('Ia_','II_')) for s in sources) if kind=='Ia_II' else any(s.startswith(('Ib_','HEEL_','TOE','TOEDF')) for s in sources)
    c.wires=[w for w in c.wires if keep(w)]
    c.records=[e for w in c.wires for e in w['edges']]
    used={e[k] for e in c.records for k in ('source','target')}
    if kind=='Ia_II':used.update(('Ia','II'))
    else:used.add('Ib')
    c.nodes={k:n for k,n in c.nodes.items() if k in used or k.startswith(('RG_','PF_','InE_','InF_'))}
    c.title='w2lvar · '+('Ia and II pathways' if kind=='Ia_II' else 'Ib, heel, and toe pathways')
    c.labels=[]
    return ORIGINAL.output(c,slug,'Same RG / Hip PF / K&A PF arrangement as the overview. Gray lines retain the RG→PF context.')

def draw_s3k():
    p,edges,d=V.load('s3k');c=ORIGINAL('s3k · RG and four PF phase channels',1250)
    V.module(c,'RG',150,['RG_E_l','RG_F_l','InE_l','InF_l'],False)
    c.group('PF phase channels',590,550,620,470)
    for name,label,x,y in [('PF_E1_l','E1',650,650),('PF_E2_l','E2',650,860),('PF_F1_l','F1',1040,650),('PF_F2_l','F2',1040,860)]:c.node(name,label,x,y,110,85,'#FFE52E')
    c.node('DRIVE','Descending drive',90,150,280,80);c.node('POSTURE','Posture command',90,340,280,80);c.node('CIN_F_r','Opposite c1-class\nCIN_F',1450,150,280,85)
    selected=[e for e in edges if e['source'] in ('DRIVE','POSTURE','CIN_F_r','RG_E_l','RG_F_l') and e['target'] in c.nodes and (not e['source'].startswith('RG_') or e['target'].startswith('PF_'))]
    count={}
    for e in selected:
        a,b=c.nodes[e['source']],c.nodes[e['target']];k=count.get(e['source'],0);count[e['source']]=k+1
        if e['source']=='CIN_F_r':points=[(1450,190),(1340,190),(1340,250),(1150,250)]
        elif e['source'].startswith('RG_'):
            left=e['source']=='RG_E_l';x=540-k*25 if left else 1260+k*25
            points=[(a['x']+55,a['y']+a['h']),(a['x']+55,400+k*30),(x,400+k*30),(x,b['y']+35),(b['x'] if left else b['x']+b['w'],b['y']+35)]
        else:
            sy=a['y']+25+k*25;ty=b['y']+20+(20 if e['source']=='POSTURE' else 0)
            if b['x']==650:points=[(a['x']+a['w'],sy),(480,sy),(480,ty),(650,ty)]
            else:points=[(a['x']+a['w'],sy),(410,sy),(410,460),(1290,460),(1290,ty),(1150,ty)]
        c.wire(points,[e],color='#888' if e['source'].startswith('RG_') else '#333')
    return c.output('s3k_architecture','Saved standing condition: actual descending, contralateral, and RG→PF connections. Internal motifs omitted.'),p

def draw_rat():
    p=ROOT/'tmp/functional_subnetwork_remote/named_conditions/rat134/Rat134_ReceptorGeometry_Standalone.asim'
    root=E.parse(p).getroot();nodes={n.findtext('ID'):n for n in root.iter('Neuron')};types={n.findtext('ID'):n for n in root.iter('SynapseType')}
    clean=lambda s:re.sub(r'\s+',' ',s.strip())
    lookup={clean(n.findtext('Name')):n.findtext('Name') for n in nodes.values()}
    edges=[]
    for e in root.iter('Connexion'):
        a,b=e.findtext('SourceID'),e.findtext('TargetID');t=types[e.findtext('SynapseTypeID')];rev=float(t.findtext('Equil'))
        edges.append({'id':e.findtext('ID'),'source':nodes[a].findtext('Name'),'target':nodes[b].findtext('Name'),'gmax_uS':float(t.findtext('SynAmp')),'reversal_mV':rev,'threshold_mV':float(t.findtext('ThreshV')),'saturation_mV':float(t.findtext('SaturateV')),'sign_category':'exc' if rev>float(nodes[b].findtext('RestingPot')) else 'inh'})
    c=ORIGINAL('Rat134 · curated receptor-geometry condition',1630)
    target_groups=[]
    for stem,label,y in [('LH_RG','RG',140),('LH_Hip PF','Hip PF',490),('LH_Knee PF','Knee PF',840),('LH_Ankle PF','Ankle PF',1190)]:
        ext=lookup[stem+' EXT'];flx=lookup[stem+' FLX'];V.module(c,label,y,[ext,flx],False);target_groups.append((ext,flx,y))
    target_set={n for a,b,y in target_groups for n in (a,b)}
    sensory=[e for e in edges if e['target'] in target_set and any(x in e['source'] for x in ('Ia','II','Ib'))]
    positions={'LH_SO- Ext Ib':(75,140),'LH_HipZ-FLX II':(75,390),'LH_IP- FLX Ia':(75,750),'LH_HipZ-Ext II':(1460,750),'Biar 1 Ia':(1460,450),'HE Ia':(1460,560)}
    aliases={'LH_SO- Ext Ib':'Ankle ext Ib (SO)','LH_HipZ-FLX II':'Hip flx II','LH_IP- FLX Ia':'Hip flx Ia (IP)','LH_HipZ-Ext II':'Hip ext II','Biar 1 Ia':'Biar 1 Ia\nprototype channel','HE Ia':'HE Ia\nprototype channel'}
    for source in sorted({e['source'] for e in sensory}):
        x,y=positions[clean(source)];c.node(source,aliases[clean(source)],x,y,275,85)
    src_count={};target_count={}
    for e in sensory:
        a,b=c.nodes[e['source']],c.nodes[e['target']];k=src_count.get(e['source'],0);src_count[e['source']]=k+1;j=target_count.get(e['target'],0);target_count[e['target']]=j+1
        sy=a['y']+18+k*18;ty=b['y']+18+j*22
        if a['x']<500:
            lane=400+({'LH_SO- Ext Ib':0,'LH_HipZ-FLX II':60,'LH_IP- FLX Ia':120}[clean(e['source'])])+k*15
            pts=[(a['x']+a['w'],sy),(lane,sy),(lane,ty),(b['x'],ty)]
        else:
            lane=1280+k*20;pts=[(a['x'],sy),(lane,sy),(lane,ty),(b['x']+b['w'],ty)]
        c.wire(pts,[e])
    for ext,flx,y in target_groups[1:]:
        feed=[e for e in edges if e['source'] in target_groups[0][:2] and e['target'] in (ext,flx)]
        lane=580 if y==490 else 1215 if y==840 else 555
        c.wire([(900,305),(900,345),(lane,345),(lane,y-35),(900,y-35),(900,y)],feed,color='#888',label='RG → PF ×2')
    return c.output('Rat134_feedback','Separate knee and ankle PF modules are retained here. These are source connections, not assumed human-model pathways.'),p

def draw_fourcell():
    p=ROOT/'Neuromechanical_Models/Mujoco_SNS_models/Hunt_fourcell_session1_left_RF_HAM_connectome.json';data=json.loads(p.read_text())
    c=ORIGINAL('RF / HAM · four-cell PF subnetwork (saved left-side projection)',1390)
    for label,x,stem in [('Hip joint PF',270,'PF_HIP'),('Knee joint PF',1010,'PF_KNEE')]:
        c.group(label,x,130,430,170)
        c.node(stem+'-E_l','ext',x+20,200,105,70,'#FFE52E');c.node(stem+'-F_l','flx',x+305,200,105,70,'#FFE52E')
    c.group('RF / HAM PF',610,540,580,440)
    c.node('PF_BI_JOINT_RECT_l','RF PF',650,620,110,85,'#FFE52E');c.node('PF_BI_JOINT_HAM_l','HAM PF',1040,620,110,85,'#FFE52E')
    c.node('PF_IN_BI_PAIR_RECT_l','RF IN',980,800,90,55,'#D9EAC1');c.node('PF_IN_BI_PAIR_HAM_l','HAM IN',760,800,90,55,'#D9EAC1')
    c.node('MN_rect_fem_l','Rectus femoris\nMN',570,1140,240,85,'#F3D5C9');c.node('MN_semimem_l','Semimembranosus\nMN',1000,1120,270,85,'#F3D5C9');c.node('MN_semiten_l','Semitendinosus\nMN',1390,1120,270,85,'#F3D5C9')
    for i,raw in enumerate(data['edges']):
        e={**raw,'id':'RFHAM:left:'+str(i),'source':raw['from'],'target':raw['to'],'gmax_uS':raw['gain'],'sign_category':raw['sign']}
        a,b=c.nodes[e['source']],c.nodes[e['target']];src=e['source'];dst=e['target']
        if src.startswith('PF_HIP') and '-E_' in src:pts=[(a['x']+52,270),(a['x']+52,370),(1095,370),(1095,620)]
        elif src.startswith('PF_HIP'):pts=[(a['x']+52,270),(a['x']+52,580),(690,580),(690,620)]
        elif src.startswith('PF_KNEE') and '-E_' in src:pts=[(a['x']+52,270),(a['x']+52,450),(725,450),(725,620)]
        elif src.startswith('PF_KNEE'):pts=[(a['x']+52,270),(a['x']+52,550),(1120,550),(1120,620)]
        elif dst=='PF_IN_BI_PAIR_RECT_l':pts=[(760,650),(875,650),(875,825),(980,825)]
        elif dst=='PF_IN_BI_PAIR_HAM_l':pts=[(1040,675),(930,675),(930,835),(850,835)]
        elif src=='PF_IN_BI_PAIR_RECT_l':pts=[(1070,825),(1165,825),(1165,680),(1150,680)]
        elif src=='PF_IN_BI_PAIR_HAM_l':pts=[(760,825),(625,825),(625,680),(650,680)]
        elif dst=='MN_rect_fem_l':pts=[(705,705),(705,1040),(690,1040),(690,1140)]
        elif dst=='MN_semimem_l':pts=[(1095,705),(1095,1040),(1135,1040),(1135,1120)]
        else:pts=[(1150,650),(1310,650),(1310,1060),(1525,1060),(1525,1120)]
        c.wire(pts,[e],label=f'{e["gmax_uS"]:.3g}')
    return c.output('RF_HAM_fourcell','Joint PF → RF/HAM pair → motor neurons. Expanded four-cell connections shown; other afferents are outside this saved projection.'),p

def hybrid_views():
    original=next(canvas for name,canvas,record in CAPTURED if name=='W2L108_feedback')
    p=ROOT/'Neuromechanical_Models/Mujoco_SNS_models/W2L_704_80_baseline/W2L_704_80_graph_ledger.json';d=json.loads(p.read_text());edge_map={e['id']:e for e in d['chemical_edges']};node_map={n['name']:n for n in d['nodes']}
    c=copy.deepcopy(original);c.title='W2L hybrid · 704 spiking and 80 graded cells (logical interfaces)';c.labels=[]
    for name,n in c.nodes.items():
        if name in node_map and node_map[name]['population']>1:n['label']+='\n'+str(node_map[name]['population'])+' cells'
    for wire in c.wires:wire['edges']=[edge_map[e['id']] for e in wire['edges']]
    c.records=[e for wire in c.wires for e in wire['edges']]
    hybrid=ORIGINAL.output(c,'W2L704_interfaces','Population readouts drive the chemical pathways; these are not soma-to-soma spike-event arrows. Internal motifs omitted.')
    p2=ROOT/'tmp/functional_subnetwork_remote/named_conditions/gait_v4/TIMING_CANDIDATE.json';cfg=json.loads(p2.read_text())
    base=next(canvas for name,canvas,record in CAPTURED if name=='syn6_feedback');c=copy.deepcopy(base);c.title='Gait2392 TimingZeroSFA · central spiking populations and PF channels'
    assert cfg['actual_spiking_cells']==256 and cfg['remaining_graded_cells']==778
    for name,n in c.nodes.items():
        if name.startswith(('RG_','PF_S')):n['label']+='\n16 cells'
    c.labels=[]
    gait=ORIGINAL.output(c,'Gait_TimingZeroSFA_interfaces','RG and S1–S6 are spiking populations; INs remain graded. Chemical pathways use spike-derived proxy voltage. Red branch remains disputed.')
    return [hybrid,gait],[p,p2]

if __name__=='__main__':
    a,p,ne=V.draw_w2l();b,vp,ve=V.draw_variant();c,sp,se=V.draw_syn6();V.key_and_findings(ne,ve,se)
    full=next(canvas for name,canvas,record in CAPTURED if name=='w2lvar_feedback')
    focused=[focus(full,'w2lvar_Ia_II','Ia_II'),focus(full,'w2lvar_Ib_contact','Ib_contact')]
    s3,source_s3=draw_s3k();rat,source_rat=draw_rat();four,source_four=draw_fourcell();hybrids,hybrid_sources=hybrid_views()
    pages=['W2L108_feedback','w2lvar_feedback','w2lvar_Ia_II','w2lvar_Ib_contact','syn6_feedback','s3k_architecture','Rat134_feedback','RF_HAM_fourcell','W2L704_interfaces','Gait_TimingZeroSFA_interfaces','synapse_key']
    combined=pymupdf.open()
    for name in pages:
        with pymupdf.open(DEST/(name+'.pdf')) as doc:combined.insert_pdf(doc)
    combined.save(DEST/'RG_PF_Schematics.pdf');combined.close()
    manifest={'pages':pages,'drawings':[a,b,c,*focused,s3,rat,four,*hybrids],'sources':{str(x):digest(x) for x in (p,vp,sp,source_s3,source_rat,source_four,*hybrid_sources)},'circuits_executed':0}
    (DEST/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({'pages':len(pages),'circuits_executed':0}))
