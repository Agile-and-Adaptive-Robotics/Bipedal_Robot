"""Representative MN drawings from saved data. Never imports circuit code."""
import copy,json,math,re
import xml.etree.ElementTree as ET
from pathlib import Path
import pymupdf
import draw_feedback_v3 as V
from audit_saved_graphs import ROOT,OUT,digest,table

DEST=OUT/'drawings/motor_examples'
DEST.mkdir(parents=True,exist_ok=True)
V.DEST=DEST
BLUE='#9AD8F1';PALE='#D8EDF7';YELLOW='#FFE52E';GRAY='#888'
PAGES=[]
class Drawing(V.Canvas):
    def node(self,name,label,x,y,w=200,h=80,color=BLUE,shape='box'):
        return super().node(name,label,x,y,w,h,color,shape)
    def edge(self,index,source,target,points,label=None):
        e=copy.deepcopy(index[source,target]);off=e.get('display_mode')=='defined_but_disabled'
        # Remove redundant collinear vertices, including zero-height routing lanes.
        compact=[]
        for p in points:
            if compact and p==compact[-1]:continue
            while len(compact)>1:
                a,b=compact[-2:];u=(b[0]-a[0],b[1]-a[1]);v=(p[0]-b[0],p[1]-b[1])
                if abs(u[0]*v[1]-u[1]*v[0])<.001 and u[0]*v[0]+u[1]*v[1]>=0:compact.pop()
                else:break
            compact.append(p)
        self.wire(compact,[e],color=GRAY if off else '#333',dashed=off,label=label)
    def output(self,slug,note=''):
        dangling=[];hidden=[];angles=[]
        for i,w in enumerate(self.wires):
            for p in (w['points'][0],w['points'][-1]):
                if not any(on_border(p,n) for n in self.nodes.values()):dangling.append((i,p))
            for a,b in zip(w['points'],w['points'][1:]):
                dx,dy=abs(a[0]-b[0]),abs(a[1]-b[1])
                if dx and dy and abs(dx-dy)>0.1:angles.append((i,a,b))
                for name,n in self.nodes.items():
                    if any(n['x']+3<a[0]+(b[0]-a[0])*t<n['x']+n['w']-3 and n['y']+3<a[1]+(b[1]-a[1])*t<n['y']+n['h']-3 for t in [k/100 for k in range(1,100)]):hidden.append((i,name))
        checks={'dangling_endpoints':dangling,'wires_hidden_by_nodes':hidden,'non_orthogonal_non_45_degree_segments':angles,'circuits_executed':0}
        (DEST/(slug+'_geometry.json')).write_text(json.dumps(checks,indent=2))
        result=super().output(slug,note)
        # MuPDF's SVG importer ignores stroke-dasharray. Expand dashes in its
        # conversion copy; keep the editable SVG's native dashed paths intact.
        svg=ET.fromstring((DEST/(slug+'.svg')).read_text(encoding='utf-8'))
        dashed=svg.findall('.//*[@stroke-dasharray]')
        if dashed:
            for path in dashed:
                path.set('d',explicit_dashes(path.attrib['d']))
                del path.attrib['stroke-dasharray']
            source=pymupdf.open(stream=ET.tostring(svg),filetype='svg')
            pdf=pymupdf.open('pdf',source.convert_to_pdf())
            pdf.save(DEST/(slug+'.pdf'))
            pdf[0].get_pixmap(matrix=pymupdf.Matrix(.85,.85)).save(DEST/(slug+'.png'))
            pdf.close();source.close()
        table(DEST/(slug+'_synapses.csv'),self.records)
        PAGES.append(slug)
        return result

def on_border(p,n):
    x,y,w,h=n['x'],n['y'],n['w'],n['h'];px,py=p
    return ((abs(px-x)<.1 or abs(px-x-w)<.1) and y-.1<=py<=y+h+.1) or ((abs(py-y)<.1 or abs(py-y-h)<.1) and x-.1<=px<=x+w+.1)

def explicit_dashes(path):
    tokens=iter(re.findall(r'[MLQA]|-?\d+(?:\.\d+)?',path));points=[];last=None
    for cmd in tokens:
        if cmd in ('M','L'):
            last=(float(next(tokens)),float(next(tokens)));points.append(last)
        elif cmd=='Q':
            control=(float(next(tokens)),float(next(tokens)));end=(float(next(tokens)),float(next(tokens)));start=last
            for k in range(1,13):
                t=k/12;points.append(tuple((1-t)**2*start[j]+2*(1-t)*t*control[j]+t*t*end[j] for j in range(2)))
            last=end
        elif cmd=='A':
            args=[float(next(tokens)) for _ in range(7)];end=tuple(args[-2:]);start=last
            center=tuple((start[j]+end[j])/2 for j in range(2));radius=V.dist(start,end)/2;angle=math.atan2(start[1]-center[1],start[0]-center[0])
            for k in range(1,19):
                a=angle+math.pi*k/18;points.append((center[0]+radius*math.cos(a),center[1]+radius*math.sin(a)))
            last=end
        else:raise ValueError(cmd)
    output=[];on=True;remaining=9.
    for start,end in zip(points,points[1:]):
        length=V.dist(start,end);used=0.
        if not length:continue
        while used<length-1e-7:
            step=min(remaining,length-used);a=tuple(start[j]+(end[j]-start[j])*used/length for j in range(2));b=tuple(start[j]+(end[j]-start[j])*(used+step)/length for j in range(2))
            if on:output.append('M '+V.pt(a)+' L '+V.pt(b))
            used+=step;remaining-=step
            if remaining<1e-7:on=not on;remaining=9. if on else 6.
    return ' '.join(output)

def index(edges):return {(e['source'],e['target']):e for e in edges}
def definitions(data):
    rows=copy.deepcopy(data['chemical_edges'])
    for i,e in enumerate(data['disabled_source_pairs']):
        e=copy.deepcopy(e)
        if not e['source'].startswith('PF_') or not e['target'].startswith('MN_'):continue
        e['id']='disabled:'+e['source']+'->'+e['target']
        e['source_defined_gmax_uS']=sum(float(t['parameters']['max_conductance']) for t in e['source_connection_definitions'])
        e['sign_category']='exc';e['display_mode']='defined_but_disabled'
        rows.append(e)
    return index(rows)

def simple_pair(native_edges):
    c=Drawing('1 · Paired MN circuit: the saved W2L108 hip subnetwork',1270);ix=index(native_edges)
    names={}
    def node(key,real,label,x,y,color=BLUE):names[key]=real;c.node(real,label,x,y,180,80,color)
    node('PE','L Hip PF ext','Hip PF\next',580,160,YELLOW);node('PF','L Hip PF flx','Hip PF\nflx',1040,160,YELLOW)
    for key,real,label,x,y in [
        ('AE','L Hip ext Ia IN','Hip ext Ia\ninput relay',180,390),('AF','L Hip flx Ia IN','Hip flx Ia\ninput relay',1440,390),
        ('IE','L Hip ext Ia','ext Ia\nreciprocal IN',580,390),('IF','L Hip flx Ia','flx Ia\nreciprocal IN',1040,390),
        ('ME','L Hip MN ext','Hip ext MN',580,670),('MF','L Hip MN flx','Hip flx MN',1040,670),
        ('BE','L Hip ext Ib','Hip ext Ib',180,670),('BF','L Hip flx Ib','Hip flx Ib',1440,670),
        ('TE','L Hip extensor II relay','Hip ext II\nrelay',180,870),('TF','L Hip flexor II relay','Hip flx II\nrelay',1440,870),
        ('RE','L Hip MN ext RE','ext RE',580,950),('RF','L Hip MN flx RE','flx RE',1040,950)]:node(key,real,label,x,y)
    def wire(a,b,p):c.edge(ix,names[a],names[b],p)
    wire('AE','IE',[(360,430),(580,430)]);wire('AF','IF',[(1440,430),(1220,430)])
    wire('PE','IE',[(670,240),(670,390)]);wire('PF','IF',[(1130,240),(1130,390)])
    wire('PE','ME',[(760,200),(830,200),(830,725),(760,725)])
    wire('PF','MF',[(1040,200),(970,200),(970,725),(1040,725)])
    wire('IE','IF',[(760,410),(1040,410)]);wire('IF','IE',[(1040,450),(760,450)])
    wire('IE','MF',[(760,460),(1000,700),(1040,700)])
    wire('IF','ME',[(1040,460),(800,700),(760,700)])
    wire('BE','ME',[(360,710),(580,710)]);wire('BF','MF',[(1440,710),(1220,710)])
    wire('TE','ME',[(360,910),(480,910),(480,690),(580,690)])
    wire('TF','MF',[(1440,910),(1320,910),(1320,690),(1220,690)])
    wire('ME','RE',[(670,750),(670,950)]);wire('MF','RF',[(1130,750),(1130,950)])
    wire('RE','ME',[(580,975),(520,975),(520,735),(580,735)])
    wire('RF','MF',[(1220,975),(1280,975),(1280,735),(1220,735)])
    wire('RE','IE',[(580,1005),(430,1005),(430,450),(580,450)])
    wire('RF','IF',[(1220,1005),(1370,1005),(1370,450),(1220,450)])
    wire('RE','RF',[(760,975),(1040,975)]);wire('RF','RE',[(1040,1005),(760,1005)])
    wire('BE','PE',[(180,710),(80,710),(80,100),(670,100),(670,160)])
    wire('IF','PF',[(1220,410),(1260,410),(1260,180),(1220,180)])
    selected={e['id'] for e in native_edges if e['source'] in c.nodes and e['target'] in c.nodes}
    assert selected=={e['id'] for e in c.records},'Missing native connections within selected subnetwork'
    c.labels += [(900,1110,'One extensor MN and one flexor MN; all saved neural connections among the displayed cells are shown.',21),
        (900,1145,'Native naming: the peripheral “Ia IN” relays feed the reciprocal-inhibition cells named “Ia”.',19),
        (900,1180,'Source detail: RE → own MN is excitatory here. RE → Ia and RE → opposite RE are inhibitory.',19)]
    return c.output('01_paired_MN','Open triangle: excitation (base on target). Filled circle: inhibition. Semicircular bridges: no junction.')

def fanout(wix):
    c=Drawing('2 · Several muscle MNs share each side of a PF pair',1050)
    c.group('K&A PF',370,130,1170,200)
    c.node('PF_KNEE-E_l','ext',410,205,240,85,YELLOW);c.node('PF_KNEE-F_l','flx',1240,205,240,85,YELLOW)
    rows=[('vas_med','Vastus medialis',90,[(460,290),(460,420),(200,420),(200,580)]),
        ('vas_int','Vastus intermedius',415,[(525,290),(525,580)]),
        ('vas_lat','Vastus lateralis',740,[(590,290),(590,460),(850,460),(850,580)]),
        ('semimem','Semimembranosus',1100,[(1300,290),(1300,420),(1210,420),(1210,580)]),
        ('semiten','Semitendinosus',1450,[(1380,290),(1380,460),(1560,460),(1560,580)])]
    for muscle,label,x,points in rows:
        mn='MN_'+muscle+'_l';pf='PF_KNEE-'+('E_l' if muscle.startswith('vas_') else 'F_l')
        c.node(mn,label+'\nMN',x,580,220,90)
        e=wix[pf,mn];c.edge(wix,pf,mn,points)
        c.labels.append((x+110,715,f'defined g = {e["source_defined_gmax_uS"]:.4g} µS',18))
    c.labels += [(900,805,'These are five separate motor pools, not one pooled “extensor” or “flexor” MN.',23),
        (900,845,'Each has its own afferents and interneurons; sharing a PF source does not merge those circuits.',21),
        (900,900,'Selected w2lvar fan-out. Additional PF inputs to the biarticular hamstrings are shown in example 3.',19),
        (900,940,'Gray dashed = source-defined connection with effective g = 0 in this saved standing condition.',19)]
    return c.output('02_shared_PF_multiple_MNs','K&A is the model’s combined knee–ankle PF pair; only the selected knee-muscle branches are expanded here.')

def biarticular(wix):
    c=Drawing('3 · A biarticular muscle receives separate inputs from two PF layers',1210)
    for label,x,stem in [('Hip PF',200,'PF_HIP'),('K&A PF',980,'PF_KNEE')]:
        c.group(label,x,130,620,200)
        c.node(stem+'-E_l','ext',x+40,210,180,85,YELLOW);c.node(stem+'-F_l','flx',x+400,210,180,85,YELLOW)
    c.node('MN_rect_fem_l','Rectus femoris\nMN',420,700,290,95)
    c.node('MN_semimem_l','Semimembranosus\nMN',1060,700,290,95)
    routes=[('PF_HIP-F_l','MN_rect_fem_l',[(690,295),(690,420),(520,420),(520,700)],(540,390)),
        ('PF_KNEE-E_l','MN_rect_fem_l',[(1110,295),(1110,550),(620,550),(620,700)],(820,525)),
        ('PF_HIP-E_l','MN_semimem_l',[(330,295),(330,480),(1160,480),(1160,700)],(780,455)),
        ('PF_KNEE-F_l','MN_semimem_l',[(1470,295),(1470,610),(1260,610),(1260,700)],(1460,660))]
    for pf,mn,pts,pos in routes:
        c.edge(wix,pf,mn,pts);c.labels.append((*pos,f'g defined = {wix[pf,mn]["source_defined_gmax_uS"]:.4g} µS',18))
    c.labels += [(565,850,'Hip flexion + knee extension',20),(1205,850,'Hip extension + knee flexion',20),
        (900,940,'Each muscle has ONE MN pool. Both PF synapses terminate on that same pool.',23),
        (900,980,'Both inputs are excitatory; separate conductances contribute to the MN membrane equation.',21),
        (900,1020,'A muscle’s action at two joints does not require a second MN or a second copy of its local reflex circuit.',20),
        (900,1100,'Gray dashed = defined source routing; every displayed PF→MN edge has effective g = 0 in the saved w2lvar condition.',19)]
    return c.output('03_biarticular_convergence','Weights belong to individual source→target synapses. Sharing an MN does not make the two PF connections one synapse.')

def local_motif(wix):
    c=Drawing('4 · One muscle’s local motif can inhibit several antagonist MNs',1390)
    nodes=[('Ia_tib_ant_l','Tibialis anterior\nIa',100,190,180,75,PALE),('II_tib_ant_l','Tibialis anterior\nII',100,490,180,75,PALE),
        ('Ib_tib_ant_l','Tibialis anterior\nIb',100,790,180,75,PALE),('IaIN_tib_ant_l','TA Ia IN',460,190,180,75,BLUE),
        ('IIIN_tib_ant_l','TA II inhibitory IN',460,410,200,75,BLUE),('IIX_tib_ant_l','TA II excitatory IN',460,610,200,75,BLUE),
        ('IBIN_tib_ant_l','TA Ib inhibitory IN',460,790,200,75,BLUE),('MN_tib_ant_l','Tibialis anterior\nMN',900,860,240,90,BLUE),
        ('RC_tib_ant_l','TA Renshaw cell',900,1090,240,75,BLUE)]
    for muscle,label,y in [('soleus','Soleus MN',190),('med_gas','Medial gastrocnemius MN',390),('lat_gas','Lateral gastrocnemius MN',590)]:nodes.append(('MN_'+muscle+'_l',label,1410,y,310,80,BLUE))
    for n in nodes:c.node(*n)
    c.edge(wix,'Ia_tib_ant_l','IaIN_tib_ant_l',[(280,227.5),(460,227.5)])
    c.edge(wix,'Ia_tib_ant_l','MN_tib_ant_l',[(190,265),(190,350),(1000,350),(1000,860)])
    c.edge(wix,'II_tib_ant_l','IIIN_tib_ant_l',[(280,510),(360,510),(360,447.5),(460,447.5)])
    c.edge(wix,'II_tib_ant_l','IIX_tib_ant_l',[(280,545),(380,545),(380,647.5),(460,647.5)])
    c.edge(wix,'IIX_tib_ant_l','MN_tib_ant_l',[(660,647.5),(780,647.5),(780,890),(900,890)])
    c.edge(wix,'Ib_tib_ant_l','IBIN_tib_ant_l',[(280,827.5),(460,827.5)])
    c.edge(wix,'IBIN_tib_ant_l','MN_tib_ant_l',[(660,827.5),(740,827.5),(740,925),(900,925)])
    c.edge(wix,'MN_tib_ant_l','RC_tib_ant_l',[(1020,950),(1020,1090)])
    c.edge(wix,'RC_tib_ant_l','MN_tib_ant_l',[(900,1110),(840,1110),(840,910),(900,910)])
    c.edge(wix,'RC_tib_ant_l','IaIN_tib_ant_l',[(900,1140),(700,1140),(700,300),(550,300),(550,265)])
    for i,(muscle,y) in enumerate([('soleus',190),('med_gas',390),('lat_gas',590)]):
        mn='MN_'+muscle+'_l';lane=1230-i*40
        sy=(215,230,247)[i]
        c.edge(wix,'IaIN_tib_ant_l',mn,[(640,sy),(lane,sy),(lane,y+25),(1410,y+25)])
        lane=1350-i*40
        sy=427.5+i*20;ty=447.5 if i==1 else y+55
        c.edge(wix,'IIIN_tib_ant_l',mn,[(660,sy),(lane,sy),(lane,ty),(1410,ty)])
    c.labels += [(900,1230,'A separate inhibitory synapse reaches each antagonist MN; the three targets are not one pooled muscle.',21),
        (900,1265,'The same TA Ia neuron also excites its own MN. II uses separate excitatory and inhibitory interneurons.',20),
        (900,1300,'PF gating, central feedback, other antagonist muscles, and the wider Renshaw network are outside this focused view.',18)]
    return c.output('04_local_reflex_fanout','TA = tibialis anterior. All lines on this page are nonzero saved w2lvar connections. No muscle or MN is duplicated.')

def nmf(syn):
    ix=index(syn['chemical_edges']);c=Drawing('5 · Active syn6 PF channels converge and fan out across muscles',1150)
    c.group('Selected channels of the six-channel PF layer',250,130,1300,200)
    for i,j in enumerate((2,3,5)):c.node(f'PF_S{j}_l',f'S{j}',390+i*400,210,220,80,YELLOW)
    for muscle,label,x in [('rect_fem','Rectus femoris\nMN',340),('vas_med','Vastus medialis\nMN',810),('soleus','Soleus\nMN',1280)]:c.node('MN_'+muscle+'_l',label,x,690,250,90)
    routes=[('PF_S2_l','MN_rect_fem_l',[(440,290),(440,690)]),
        ('PF_S2_l','MN_vas_med_l',[(550,290),(550,420),(890,420),(890,690)]),
        ('PF_S3_l','MN_rect_fem_l',[(840,290),(840,490),(520,490),(520,690)]),
        ('PF_S3_l','MN_soleus_l',[(950,290),(950,450),(1360,450),(1360,690)]),
        ('PF_S5_l','MN_rect_fem_l',[(1240,290),(1240,550),(570,550),(570,690)]),
        ('PF_S5_l','MN_soleus_l',[(1350,290),(1350,380),(1450,380),(1450,690)])]
    for a,b,pts in routes:c.edge(ix,a,b,pts)
    selected=[e for e in syn['chemical_edges'] if e['target'] in c.nodes and e['source'] in c.nodes]
    assert {e['id'] for e in selected}=={e['id'] for e in c.records},'Selected PF-channel outputs must be complete'
    c.labels += [(900,845,'All three PF inputs to rectus femoris are shown. The other two MNs also receive inputs from channels outside this view.',19),
        (900,890,'Rectus femoris integrates three separate PF synapses; S2 also branches to vastus medialis.',22),
        (900,935,'S3 and S5 also branch to soleus. These channels encode muscle synergies, not named hip/knee half-centers.',20),
        (900,990,'S2 → rectus femoris: g = 0.007152 µS; S2 → vastus medialis: g = 0.01747 µS.',20),
        (900,1035,'Each outgoing synapse has its own weight. Other channel inputs and local reflex circuits are omitted here.',19)]
    return c.output('05_active_synergy_convergence','All displayed connections are excitatory and nonzero. Exact conductances and thresholds accompany the editable drawing.')

def verify(sources,records):
    lookup={}
    for source,rows in sources:
        for e in rows:lookup[e['id']]=(source,e)
    reports=[]
    for slug in PAGES:
        d=json.loads((DEST/(slug+'_connections.json')).read_text())
        for e in d['source_connections']:
            source,original=lookup[e['id']]
            assert e==original,(slug,e['id'])
        g=json.loads((DEST/(slug+'_geometry.json')).read_text())
        assert not any(g[k] for k in ('dangling_endpoints','wires_hidden_by_nodes','non_orthogonal_non_45_degree_segments')),(slug,g)
        reports.append({'page':slug,'saved_source_records':len(d['source_connections']),'geometry_pass':True,'exact_source_record_parity':True})
    report={'pages':reports,'sources':{str(p):digest(p) for p,rows in sources},'circuits_executed':0,'scope':'Selected displayed saved connections; no behavioral or complete-model claim.'}
    (DEST/'VERIFICATION.json').write_text(json.dumps(report,indent=2))
    return report

if __name__=='__main__':
    np,ne=V.native();wp,we,wd=V.load('w2lvar');sp,se,sd=V.load('syn6');wix=definitions(wd)
    records=[simple_pair(ne),fanout(wix),biarticular(wix),local_motif(wix),nmf(sd)]
    report=verify([(np,ne),(wp,list(wix.values())),(sp,se)],records)
    packet=pymupdf.open()
    for slug in PAGES:
        with pymupdf.open(DEST/(slug+'.pdf')) as page:packet.insert_pdf(page)
    titles=['Paired W2L108 MN circuit','Shared PF: several agonist and antagonist MNs','Biarticular: Hip and K&A PF inputs to one MN','Local reflex fan-out to several antagonist MNs','Active syn6 convergence and shared PF outputs']
    packet.set_toc([[1,title,i+1] for i,title in enumerate(titles)])
    packet.save(DEST/'Representative_MN_Layers.pdf');packet.close()
    combined=pymupdf.open(OUT/'drawings/v4/RG_PF_Schematics.pdf')
    previous_pages=len(combined)
    with pymupdf.open(DEST/'Representative_MN_Layers.pdf') as motor:combined.insert_pdf(motor)
    old=json.loads((OUT/'drawings/v4/MANIFEST.json').read_text())['pages']
    combined.set_toc([[1,'RG and PF schematics',1]]+[[2,name.replace('_',' '),i+1] for i,name in enumerate(old)]+[[1,'Representative MN layers',previous_pages+1]]+[[2,title,previous_pages+i+1] for i,title in enumerate(titles)])
    combined.save(DEST/'RG_PF_and_Representative_MN_Layers.pdf');combined.close()
    print(json.dumps({'motor_pages':len(PAGES),'combined_pages':11+len(PAGES),'source_connection_records':sum(p['saved_source_records'] for p in report['pages']),'circuits_executed':0}))
