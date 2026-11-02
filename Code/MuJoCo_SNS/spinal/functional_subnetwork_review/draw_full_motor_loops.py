"""Draw saved neural connections and physical maps; never execute a circuit."""
import copy,json,math
from pathlib import Path
import xml.etree.ElementTree as E
import pymupdf
import draw_motor_examples as M
import draw_feedback_v3 as V
from audit_saved_graphs import OUT,table,digest

DEST=OUT/'drawings/full_motor_loops';DEST.mkdir(parents=True,exist_ok=True)
M.DEST=DEST;V.DEST=DEST
WIDTH=2400;PHYSICAL='#9E5F45';SALMON='#FA8E80'

class FullDrawing(M.Drawing):
    def __init__(self,title):
        super().__init__(title,1900);self.physical_records=[];self.icon_nodes={};self.aliases={}
    def physical(self,source,target,points,policy,units):
        record={'id':'map:'+source+'->'+target,'source':source,'target':target,'kind':'physical_signal','policy':policy,'units':units}
        self.wire(points,[record],sign='physical',color=PHYSICAL)
        self.physical_records.append(record)
    def output(self,slug,note=''):
        result=super().output(slug,'')
        path=DEST/(slug+'.svg');root=E.fromstring(path.read_text(encoding='utf-8'));ns='{http://www.w3.org/2000/svg}'
        root.set('width',str(WIDTH));root.set('viewBox',f'0 0 {WIDTH} {self.h}')
        root.find(ns+'text').set('x',str(WIDTH/2))
        # Physical signals are directed arrows, not excitatory synapses.
        for elem in list(root):
            if elem.tag==ns+'polygon' and elem.get('stroke')==PHYSICAL:root.remove(elem)
        for wire in self.wires:
            if wire['sign']!='physical':continue
            end=wire['points'][-1];before=wire['points'][-2];back=V.toward(end,before,11)
            length=V.dist(before,end);ux=(end[0]-before[0])/length;uy=(end[1]-before[1])/length
            a=(back[0]-uy*5,back[1]+ux*5);b=(back[0]+uy*5,back[1]-ux*5)
            root.append(E.Element(ns+'polyline',{'points':f'{V.pt(a)} {V.pt(end)} {V.pt(b)}','fill':'none','stroke':PHYSICAL,'stroke-width':'1.8','stroke-linejoin':'round'}))
        # Vector muscle, spindle, and tendon-force symbols; no model changes.
        for name,kind in self.icon_nodes.items():
            n=self.nodes[name];x,y,w,h=[n[k] for k in ('x','y','w','h')];cx=x+w/2;cy=y+h/2
            rect=next(e for e in root if e.tag==ns+'rect' and e.get('x')==str(x) and e.get('y')==str(y) and e.get('width')==str(w))
            if kind=='muscle':
                root.remove(rect)
                root.append(E.Element(ns+'path',{'d':f'M {cx},{y} C {x+w+w/6},{y+h*.17} {x+w+w/6},{y+h*.83} {cx},{y+h} C {x-w/6},{y+h*.83} {x-w/6},{y+h*.17} {cx},{y} Z','fill':'#F0C7AE','stroke':'#654C3D','stroke-width':'2'}))
                for dx in (-20,20):root.append(E.Element(ns+'path',{'d':f'M {cx},{y+10} Q {cx+dx},{cy} {cx},{y+h-10}','fill':'none','stroke':'#B78A72','stroke-width':'1'}))
                label=n.get('muscle_label','Muscle');root.extend(E.fromstring('<g xmlns="http://www.w3.org/2000/svg">'+V.text(cx,cy+7,label,21)+'</g>'))
            elif kind=='spindle':
                coil=[(x+20,cy)]
                for k in range(11):coil.append((x+30+k*(w-60)/10,cy+(-10 if k%2==0 else 10)))
                coil.append((x+w-20,cy))
                root.append(E.Element(ns+'polyline',{'points':' '.join(V.pt(p) for p in coil),'fill':'none','stroke':'#775A43','stroke-width':'1.6'}))
                root.extend(E.fromstring('<g xmlns="http://www.w3.org/2000/svg">'+V.text(cx,y+20,'Spindle',17)+V.text(cx,y+h-8,'l, dl/dt',17)+'</g>'))
            elif kind=='tendon':
                root.remove(rect)
                root.append(E.Element(ns+'path',{'d':f'M {cx},{y} L {x+w},{cy} L {cx},{y+h} L {x},{cy} Z','fill':'#F8E2CF','stroke':'#775A43','stroke-width':'1.6'}))
                root.extend(E.fromstring('<g xmlns="http://www.w3.org/2000/svg">'+V.text(cx,cy+6,'Force / Ib',17)+'</g>'))
            elif kind=='motor_bundle':
                root.extend(E.fromstring('<g xmlns="http://www.w3.org/2000/svg">'+V.text(cx,y+28,'4 separate MN pools',19)+'</g>'))
                for label,xx,yy in [('SM',x+15,y+50),('ST',x+150,y+50),('BFLH',x+15,y+120),('BFSH',x+150,y+120)]:
                    root.append(E.Element(ns+'rect',{'x':str(xx),'y':str(yy),'width':'115','height':'45','rx':'10','fill':M.BLUE,'stroke':'#333','stroke-width':'1.3'}))
                    root.extend(E.fromstring('<g xmlns="http://www.w3.org/2000/svg">'+V.text(xx+57.5,yy+29,label,19)+'</g>'))
        svg=E.tostring(root,encoding='unicode');path.write_text(svg,encoding='utf-8')
        # Preserve visible dashes in the PDF conversion as well as editable SVG.
        for elem in root.findall('.//*[@stroke-dasharray]'):
            elem.set('d',M.explicit_dashes(elem.get('d')));del elem.attrib['stroke-dasharray']
        src=pymupdf.open(stream=E.tostring(root),filetype='svg');doc=pymupdf.open('pdf',src.convert_to_pdf())
        doc.save(DEST/(slug+'.pdf'));doc[0].get_pixmap(matrix=pymupdf.Matrix(.65,.65)).save(DEST/(slug+'.png'));doc.close();src.close()
        neural=[r for r in self.records if r.get('kind')!='physical_signal']
        (DEST/(slug+'_connections.json')).write_text(json.dumps({'title':self.title,'neural_connections':neural,'physical_signal_paths':self.physical_records,'display_bundles':self.aliases,'circuits_executed':0},indent=2))
        table(DEST/(slug+'_synapses.csv'),neural);table(DEST/(slug+'_physical_maps.csv'),self.physical_records)
        return result

def make_unit(c,muscle,short,right=False):
    def point(x,y):return (WIDTH-50-x if right else 50+x,y)
    def points(p):return [point(x,y) for x,y in p]
    def node(key,label,x,y,w,h,color=M.BLUE,shape='box'):
        xx=WIDTH-50-x-w if right else 50+x
        name=(key+'_'+muscle+'_l') if key in ('Ia','II','Ib','IaIN','IIX','IIIN','IBIN','MN','RC') else key+':'+muscle
        c.node(name,label,xx,y,w,h,color,shape);return name
    names={}
    for key,label,x,y,w,h,color,shape in [
        ('Ia','Ia',400,450,160,70,M.PALE,'box'),('II','II',400,600,160,70,M.PALE,'box'),('Ib','Ib',400,1450,160,70,M.PALE,'box'),
        ('IaIN','Ia IN',650,450,200,70,M.BLUE,'box'),('IIX','II excitatory IN',650,600,200,70,M.BLUE,'box'),('IIIN','II inhibitory IN',650,750,200,70,M.BLUE,'box'),
        ('IBIN','Ib inhibitory IN',650,1450,200,70,M.BLUE,'box'),('MN','MN',650,1150,200,100,M.BLUE,'box'),('RC','Renshaw',650,1320,200,80,M.BLUE,'box'),
        ('muscle','',70,1130,150,140,'#F0C7AE','box'),('spindle','',70,450,150,90,'#FFF1D2','box'),('tendon','',70,1450,150,70,'#F8E2CF','box'),
        ('AD_Ia','AD',280,450,70,70,SALMON,'circle'),('AD_II','AD',280,600,70,70,SALMON,'circle'),('AD_Ib','AD',280,1450,70,70,SALMON,'circle'),('AD_motor','AD',315,1160,80,80,SALMON,'circle')]:
        names[key]=node(key,label,x,y,w,h,color,shape)
    for key in ('muscle','spindle','tendon'):c.icon_nodes[names[key]]=key
    c.nodes[names['muscle']]['muscle_label']=short
    full_names={'vas_med':'Vastus medialis','vas_lat':'Vastus lateralis','rect_fem':'Rectus femoris','semimem':'Semimembranosus'}
    c.labels.append((*point(300,365),short+' · '+full_names[muscle],25))
    physical=[('MN','AD_motor',[(650,1200),(395,1200)],'clip(MN voltage / 5 mV,0,1)','mV to activation'),
        ('AD_motor','muscle',[(315,1200),(220,1200)],'native activation and contraction dynamics','activation'),
        ('muscle','spindle',[(145,1130),(145,540)],'muscle-tendon length and length rate supplied to spindle-proxy encoders','m; m/s'),
        ('spindle','AD_Ia',[(220,485),(280,485)],'Ia length-rate encoder; saved AFF/MOD/G gating','m/s'),
        ('AD_Ia','Ia',[(350,485),(400,485)],'saved Ia input injection; weight 1','nA'),
        ('spindle','AD_II',[(220,515),(250,515),(250,635),(280,635)],'II length encoder relative to initial standing length; saved gating','m'),
        ('AD_II','II',[(350,635),(400,635)],'saved II input injection; weight 1','nA'),
        ('muscle','tendon',[(145,1270),(145,1450)],'abs(native muscle force)','N'),
        ('tendon','AD_Ib',[(220,1485),(280,1485)],'Ib force/Fmax encoder; saved AFF/MOD/G gating','N'),
        ('AD_Ib','Ib',[(350,1485),(400,1485)],'saved Ib input injection; weight 1','nA')]
    for a,b,p,policy,units in physical:c.physical(names[a],names[b],points(p),policy,units)
    local=[('Ia','IaIN',[(560,485),(650,485)]),('Ia','MN',[(480,520),(480,560),(610,560),(610,1180),(650,1180)]),
        ('II','IIX',[(560,635),(650,635)]),('II','IIIN',[(560,655),(590,655),(590,785),(650,785)]),
        ('IIX','MN',[(850,635),(890,635),(890,1200),(850,1200)]),('Ib','IBIN',[(560,1485),(650,1485)]),
        ('IBIN','MN',[(750,1450),(750,1420),(590,1420),(590,1210),(650,1210)]),
        ('MN','RC',[(750,1250),(750,1320)]),('RC','MN',[(650,1360),(600,1360),(600,1230),(650,1230)]),
        ('RC','IaIN',[(850,1360),(950,1360),(950,560),(750,560),(750,520)])]
    return names,points,local

def figure(ix,kind):
    bi=kind=='biarticular';muscles=['rect_fem','semimem'] if bi else ['vas_med','vas_lat'];shorts=['RF','SM'] if bi else ['VM','VL']
    c=FullDrawing('Biarticular muscles · PF inputs, local reflexes, and physical maps' if bi else 'Shared PF source · two muscle-specific reflex and physical loops')
    units=[make_unit(c,m,s,bool(i)) for i,(m,s) in enumerate(zip(muscles,shorts))]
    for names,transform,local in units:
        for a,b,p in local:c.edge(ix,names[a],names[b],transform(p))
    left,right=units[0][0],units[1][0]
    # Renshaw cross-inhibition is present even for the two agonists.
    for a,b,y in [(left['RC'],right['RC'],1340),(right['RC'],left['RC'],1380)]:
        start=(900,y) if a==left['RC'] else (1500,y);end=(1500,y) if a==left['RC'] else (900,y)
        c.edge(ix,a,b,[start,end])
    if bi:
        for key,y1,y2 in [('IaIN',470,495),('IBIN',1470,1500)]:
            c.edge(ix,left[key],right[key],[(900,y1),(1500,y1)])
            c.edge(ix,right[key],left[key],[(1500,y2),(900,y2)])
        c.edge(ix,left['IaIN'],right['MN'],[(900,505),(960,505),(1480,1025),(1580,1025),(1580,1150)])
        c.edge(ix,right['IaIN'],left['MN'],[(1500,505),(1440,505),(920,1025),(820,1025),(820,1150)])
        c.edge(ix,left['IIIN'],right['MN'],[(900,785),(1050,785),(1300,1035),(1620,1035),(1620,1150)])
        c.edge(ix,right['IIIN'],left['MN'],[(1500,785),(1350,785),(1100,1035),(780,1035),(780,1150)])
    else:
        bundle='Knee-flexor MN destinations';targets=['MN_'+m+'_l' for m in ('semimem','semiten','bifemlh','bifemsh')]
        c.node(bundle,'',1060,1040,280,190,'#F8FBFD');c.icon_nodes[bundle]='motor_bundle';c.aliases[bundle]=targets
        for source,pts in [
            (left['IaIN'],[(900,505),(1080,505),(1080,1000),(1150,1000),(1150,1040)]),
            (right['IaIN'],[(1500,505),(1320,505),(1320,1000),(1250,1000),(1250,1040)]),
            (left['IIIN'],[(900,785),(980,785),(980,1150),(1060,1150)]),
            (right['IIIN'],[(1500,785),(1420,785),(1420,1150),(1340,1150)])]:
            rows=[ix[source,t] for t in targets]
            assert len({(e['gmax_uS'],e['reversal_mV'],e['threshold_mV'],e['saturation_mV']) for e in rows})==1
            c.wire(pts,rows)
        c.labels.append((1200,1285,'Each line = ×4 synapses',18))
    pf=[]
    if bi:
        for title,x,stem in [('Hip PF',480,'PF_HIP'),('K&A PF',1260,'PF_KNEE')]:
            c.group(title,x,100,660,170)
            for suffix,xx in [('E',x+70),('F',x+400)]:
                name=stem+'-'+suffix+'_l';c.node(name,'ext' if suffix=='E' else 'flx',xx,165,190,85,M.YELLOW);pf.append(name)
        inputs=[('PF_HIP-F_l',left['MN'],[(975,250),(975,290),(1020,290),(1020,1168),(900,1168)]),
            ('PF_KNEE-E_l',left['MN'],[(1385,250),(1385,315),(1060,315),(1060,1240),(900,1240)]),
            ('PF_HIP-E_l',right['MN'],[(645,250),(645,285),(1320,285),(1320,1168),(1500,1168)]),
            ('PF_KNEE-F_l',right['MN'],[(1755,250),(1755,325),(1360,325),(1360,1240),(1500,1240)])]
        gate='PF_KNEE-F_l';starts=[(1700,250),(1810,250)]
    else:
        c.group('K&A PF',810,100,780,170)
        for suffix,x in [('E',880),('F',1330)]:
            name='PF_KNEE-'+suffix+'_l';c.node(name,'ext' if suffix=='E' else 'flx',x,165,190,85,M.YELLOW);pf.append(name)
        inputs=[('PF_KNEE-E_l',left['MN'],[(925,250),(925,300),(1020,300),(1020,1170),(900,1170)]),
            ('PF_KNEE-E_l',right['MN'],[(1025,250),(1025,285),(1360,285),(1360,1170),(1500,1170)])]
        gate='PF_KNEE-F_l';starts=[(1375,250),(1475,250)]
    for a,b,p in inputs:c.edge(ix,a,b,p)
    if not bi:
        rows=[ix['PF_KNEE-F_l',t] for t in targets]
        c.wire([(1425,250),(1425,405),(1200,405),(1200,1040)],rows,color=M.GRAY,dashed=True)
    for i,(unit,start) in enumerate(zip(units,starts)):
        target=unit[0]['IaIN'];lane=1100 if i==0 else 1300;y=340+i*35;end=(900,475) if i==0 else (1500,475)
        c.edge(ix,gate,target,[start,(start[0],y),(lane,y),(lane,475),end])
    shown={n for n in c.nodes if n.startswith(('PF_','MN_','RC_','Ia_','II_','Ib_','IaIN_','IIX_','IIIN_','IBIN_'))}
    for members in c.aliases.values():shown.update(members)
    central_feedback=[e for e in ix.values() if e['source'] in shown and e['target'] in shown and e['source'].startswith(('Ia_','II_','Ib_')) and e['target'].startswith('PF_')]
    central_ids={e['id'] for e in central_feedback}
    expected={e['id'] for e in ix.values() if e['source'] in shown and e['target'] in shown and e['id'] not in central_ids}
    actual={e['id'] for e in c.records if e.get('kind')!='physical_signal'}
    assert expected==actual,{'missing':expected-actual,'extra':actual-expected}
    c.labels += [(1200,1620,'AD = map / adapter. Brown arrows: physical signals. Open triangles: excitation. Filled circles: inhibition.',23),
        (1200,1660,'Spindle and tendon symbols show encoded signals; muscle→spindle→AD→Ia/II and muscle→force sensor→AD→Ib are explicit.',21),
        (1200,1700,'MN→AD→muscle carries activation. Each muscle retains its own MN, Ia, II, Ib, interneurons, and Renshaw cell.',22),
        (1200,1750,'Gray dashed PF→MN routes have effective g = 0. Central afferent→RG/PF projections are detailed in the preceding RG/PF sheets.',20)]
    c.labels.append((1200,1795,'RF: hip flexion + knee extension. SM: hip extension + knee flexion. Both PF inputs terminate on the same muscle MN.' if bi else 'VM and VL share PF drive, but retain separate loops. The four antagonist destinations are grouped; every line to that group represents ×4 synapses.',20))
    if not bi:c.labels.append((1200,1835,'SM = semimembranosus; ST = semitendinosus; BFLH/BFSH = biceps femoris long/short head. Their own local loops are outside this crop.',19))
    slug='biarticular_full_loops' if bi else 'multiple_muscles_full_loops';c.output(slug)
    return slug,c

def main():
    source,edges,data=V.load('w2lvar');ix=M.definitions(data)
    results=[figure(ix,'multiple'),figure(ix,'biarticular')]
    reports=[]
    for slug,c in results:
        for e in c.records:
            if e.get('kind')!='physical_signal':assert e==ix[e['source'],e['target']]
        for e in c.physical_records:
            if e['target'].startswith(('Ia_','II_','Ib_')):
                candidates=[a for a in data['input_adapters'] if a['input']==e['target'] and a['target']==e['target']]
                assert len(candidates)==1 and candidates[0]['weight']==1
        checks=json.loads((DEST/(slug+'_geometry.json')).read_text())
        assert not any(checks[k] for k in ('dangling_endpoints','wires_hidden_by_nodes','non_orthogonal_non_45_degree_segments')),(slug,checks)
        reports.append({'drawing':slug,'neural_records':len(c.records)-len(c.physical_records),'physical_arrows':len(c.physical_records),'all_local_and_descending_connections_between_displayed_neurons':True,'physical_input_ports_verified':True,'geometry_pass':True,'central_afferent_to_RG_PF_paths':'See preceding RG/PF sheets; this view expands local motor circuits.'})
    (DEST/'VERIFICATION.json').write_text(json.dumps({'drawings':reports,'source':str(source),'source_sha256':digest(source),'recorded_actuation':data['actuation'],'recorded_encoder_policy':data['encoder_units_and_signs'],'circuits_executed':0},indent=2))
    expanded=pymupdf.open()
    for slug,c in results:
        with pymupdf.open(DEST/(slug+'.pdf')) as doc:expanded.insert_pdf(doc)
    expanded.set_toc([[1,'Multiple muscles: full local loops and physical maps',1],[1,'Biarticular pair: full coupled loops and physical maps',2]])
    expanded.save(DEST/'Expanded_MN_Feedback_and_Muscle_Maps.pdf');expanded.close()
    old=OUT/'drawings/motor_examples';files=[old/'01_paired_MN.pdf',DEST/'multiple_muscles_full_loops.pdf',DEST/'biarticular_full_loops.pdf',old/'04_local_reflex_fanout.pdf',old/'05_active_synergy_convergence.pdf']
    packet=pymupdf.open()
    for path in files:
        with pymupdf.open(path) as doc:packet.insert_pdf(doc)
    titles=['Paired W2L108 MN circuit','Multiple muscles: full loops and muscle maps','Biarticular pair: full loops and muscle maps','Local reflex fan-out','Active syn6 convergence']
    packet.set_toc([[1,t,i+1] for i,t in enumerate(titles)]);packet.save(DEST/'Representative_MN_Layers_Expanded.pdf');packet.close()
    combined=pymupdf.open(OUT/'drawings/v4/RG_PF_Schematics.pdf')
    with pymupdf.open(DEST/'Representative_MN_Layers_Expanded.pdf') as doc:combined.insert_pdf(doc)
    combined.set_toc([[1,'RG and PF schematics',1],[1,'Representative MN layers',12]]+[[2,t,i+12] for i,t in enumerate(titles)])
    combined.save(DEST/'RG_PF_MN_Full_Feedback.pdf');combined.close()
    print(json.dumps({'drawings':reports,'combined_pages':16,'circuits_executed':0}))

if __name__=='__main__':main()
