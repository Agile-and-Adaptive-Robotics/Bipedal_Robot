"""Literature-style RG/PF feedback pages from the unchanged saved W2L108 XML."""
import html,json,math
from pathlib import Path
import xml.etree.ElementTree as ET
import pymupdf
from audit_saved_graphs import ROOT,OUT,table,digest

DEST=OUT/'drawings';DEST.mkdir(exist_ok=True)
SOURCE=ROOT/'tmp/functional_subnetwork_remote/W2L108_verified/editable/W2L_VirtualWalker_Newer_Standalone.asim'
r=ET.parse(SOURCE).getroot();nodes={n.findtext('ID'):n for n in r.iter('Neuron')};types={t.findtext('ID'):t for t in r.iter('SynapseType')}
edges=[]
for e in r.iter('Connexion'):
    a,b=e.findtext('SourceID'),e.findtext('TargetID');t=types[e.findtext('SynapseTypeID')]
    edges.append({'id':e.findtext('ID'),'source_id':a,'target_id':b,'source':nodes[a].findtext('Name'),'target':nodes[b].findtext('Name'),'type_id':e.findtext('SynapseTypeID'),'type_name':t.findtext('Name'),'gmax_uS':float(t.findtext('SynAmp')),'Erev_mV':float(t.findtext('Equil')),'threshold_mV':float(t.findtext('ThreshV')),'saturation_mV':float(t.findtext('SaturateV')),'target_rest_mV':float(nodes[b].findtext('RestingPot'))})

def text(x,y,s,size=17,anchor='middle',fill='#111'):
    lines=s.split('\n');return ''.join(f'<text x="{x}" y="{y+i*size*1.25}" font-family="Arial" font-size="{size}" text-anchor="{anchor}" fill="{fill}">{html.escape(line)}</text>' for i,line in enumerate(lines))
def rect(x,y,w,h,fill):return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{fill}" stroke="#222" stroke-width="1.7"/>'
def start(w,h):return f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><defs><marker id="exc" viewBox="0 0 12 12" refX="11" refY="6" markerWidth="8" markerHeight="8" orient="auto-start-reverse"><path d="M1,1 L11,6 L1,11 Z" fill="white" stroke="#333" stroke-width="1.2"/></marker><marker id="inh" viewBox="0 0 12 12" refX="11" refY="6" markerWidth="8" markerHeight="8" orient="auto"><circle cx="6" cy="6" r="4.5" fill="#111"/></marker></defs>'
def save(slug,svg):
    p=DEST/(slug+'.svg');p.write_text(svg,encoding='utf-8')
    doc=pymupdf.open(stream=svg.encode(),filetype='svg');pdf=pymupdf.open('pdf',doc.convert_to_pdf());pdf.save(DEST/(slug+'.pdf'));pdf[0].get_pixmap(matrix=pymupdf.Matrix(1.5,1.5)).save(DEST/(slug+'.png'));pdf.close();doc.close()

for side in ('L','R'):
    contra='R' if side=='L' else 'L'
    targets=[side+' RG ext',side+' RG flx',side+' Hip PF ext',side+' Hip PF flx',side+' Knee PF ext',side+' Knee PF flx']
    in_targets=[n+' IN' for n in targets]
    selected=[e for e in edges if (e['target'] in targets+in_targets) and (any(k in e['source'] for k in ('Ia','Ib','II','heel','toe')) or e['source']==contra+' V3')]
    # The target ports below are actual source neuron names. Internal RG/PF cells
    # stay inside their subsystem, as in Ben's supplied reference image.
    modules={side+' RG':190,side+' Hip PF':550,side+' Knee PF':910}
    sensors={side+' Hip ext Ib':(35,210),side+' Knee ext Ib':(35,550),side+' heel contact':(35,910),side+' Hip flx Ia':(1130,125),side+' Hip flx II':(1130,345),side+' Knee flx Ia':(1130,565),side+' Knee flx II':(1130,785),contra+' heel contact':(1130,1005),contra+' V3':(1130,1225)}
    svg=start(1400,1530)+text(700,38,f'W2L108 · {"left" if side=="L" else "right"} central sensory feedback',25)
    svg+=text(700,68,'Half-centres and interneurons shown explicitly; one source neuron per sensory box',18)
    # Route each actual feedback edge on its own track. Connections are drawn
    # before symbols; crossings have no junction dots and never imply a synapse.
    rows=[];counts={};tips=[]
    for i,e in enumerate(selected):
        code=f'{side}{i+1:02d}';e=dict(e,figure_key=code);rows.append(e)
        x,y=sensors[e['source']];left=x<400;startx=x+235 if left else x;starty=y+32
        target_name=e['target'].removesuffix(' IN');stem=target_name.rsplit(' ',1)[0];my=modules[stem];ext=target_name.endswith('ext')
        targetx=500 if ext else 900;targety=my+53
        key=(stem,ext);slot=counts.get(key,0);counts[key]=slot+1;targety+=slot*7-10
        lane=(325+i*7) if left else (1090-i*7)
        marker='exc' if e['Erev_mV']>e['target_rest_mV'] else 'inh'
        if e['target'].endswith(' IN'):
            targetx=624 if ext else 744;targety=my+120
            points=f'{startx},{starty} {lane},{starty} {lane},{my+105} {targetx},{my+105} {targetx},{targety}'
        else:points=f'{startx},{starty} {lane},{starty} {lane},{targety} {targetx},{targety}'
        svg+=f'<g id="{e["id"]}"><title>{html.escape(e["source"]+" → "+e["target"]+" | "+str(e["gmax_uS"])+" uS | "+e["type_name"])}</title><polyline points="{points}" fill="none" stroke="#333" stroke-width="1.5" stroke-linejoin="round"/></g>'
        direction=1 if ext else -1
        if e['target'].endswith(' IN'):
            tips.append(f'<polygon points="{targetx},{targety-12} {targetx-6},{targety} {targetx+6},{targety}" fill="white" stroke="#222" stroke-width="1.5"/>' if marker=='exc' else f'<circle cx="{targetx}" cy="{targety-5}" r="5" fill="#111"/>')
        else:tips.append(f'<polygon points="{targetx-12*direction},{targety} {targetx},{targety-6} {targetx},{targety+6}" fill="white" stroke="#222" stroke-width="1.5"/>' if marker=='exc' else f'<circle cx="{targetx-5*direction}" cy="{targety}" r="5" fill="#111"/>')
        svg+=text(lane-4 if left else lane+4,(starty+targety)/2,code,14,'end' if left else 'start', '#555')
    # Exact RG→PF output edges, with separate extensor/flexor lanes.
    feed=[e for e in edges if e['source'] in (side+' RG ext',side+' RG flx') and e['target'] in targets and 'PF' in e['target']]
    for e in feed:
        ext=e['source'].endswith('ext');x=545 if ext else 855;my=modules[e['target'].rsplit(' ',1)[0]]
        sx=x;sy=modules[side+' RG']+95;targety=my+40;lane=460 if ext else 940
        svg+=f'<polyline points="{sx},{sy} {lane},{sy+90} {lane},{targety-25} {x},{targety-25} {x},{targety}" fill="none" stroke="#777" stroke-width="1.7"/>'
        tips.append(f'<polygon points="{x},{targety-12} {x-6},{targety} {x+6},{targety}" fill="white" stroke="#777" stroke-width="1.5"/>')
    for stem,y in modules.items():
        svg+=rect(500,y,400,180,'none')+rect(590,y,220,32,'#91E714')+rect(500,y+40,90,55,'#FFE52E')+rect(810,y+40,90,55,'#FFE52E')
        svg+=text(545,y+63,'ext\nHC',16)+text(855,y+63,'flx\nHC',16)
        title='RG' if stem.endswith(' RG') else 'Hip PF' if 'Hip' in stem else 'Knee / ankle\nPF'
        svg+=text(700,y+22,title.replace('\n',' '),18)
        svg+=rect(600,y+120,80,45,'#D3EABD')+rect(720,y+120,80,45,'#D3EABD')
        svg+=text(640,y+139,'ext IN',15)+text(760,y+139,'flx IN',15)
        # These four internal edges are verified source edges, not inferred wiring.
        internal=[e for e in edges if e['source'].startswith(stem+' ') and e['target'].startswith(stem+' ') and e['source'] in targets+in_targets and e['target'] in targets+in_targets]
        for e in internal:
            if e['source']==stem+' ext' and e['target']==stem+' ext IN':
                pts=f'590,{y+65} 640,{y+65} 640,{y+120}';tip=f'<polygon points="640,{y+108} 634,{y+120} 646,{y+120}" fill="white" stroke="#333"/>'
            elif e['source']==stem+' flx' and e['target']==stem+' flx IN':
                pts=f'810,{y+65} 760,{y+65} 760,{y+120}';tip=f'<polygon points="760,{y+108} 754,{y+120} 766,{y+120}" fill="white" stroke="#333"/>'
            elif e['source']==stem+' ext IN' and e['target']==stem+' flx':
                pts=f'680,{y+145} 690,{y+145} 690,{y+175} 920,{y+175} 920,{y+84} 900,{y+84}';tip=f'<circle cx="905" cy="{y+84}" r="5" fill="black"/>'
            elif e['source']==stem+' flx IN' and e['target']==stem+' ext':
                pts=f'720,{y+155} 700,{y+155} 700,{y+185} 480,{y+185} 480,{y+84} 500,{y+84}';tip=f'<circle cx="495" cy="{y+84}" r="5" fill="black"/>'
            else:continue
            svg+=f'<g id="{e["id"]}"><polyline points="{pts}" fill="none" stroke="#333" stroke-width="1.5"/>{tip}</g>'
        feed+=internal
    for name,(x,y) in sensors.items():
        svg+=rect(x,y,235,64,'#FFECA8')
        s=name.replace(side+' ','',1) if name.startswith(side+' ') else name+'\n(opposite side)'
        svg+=text(x+117.5,y+27,s,17)
    svg+=''.join(tips)
    svg+=rect(35,1170,235,95,'#EEE')+text(152,1198,'Toe contact\nNo direct RG/PF edge\nin this condition',16)
    svg+=text(700,1350,'Triangle base meets its receiving object; triangle point lies on the incoming line.',17)
    svg+=text(700,1380,'Filled circle: inhibition at target rest. Crossings without dots are not junctions.',17)
    svg+=text(700,1410,'Ia / II / Ib sensory boxes each represent one source afferent neuron, not a muscle bundle.',17)
    svg+=text(700,1440,'All direct afferent targets here are HCs; opposite-side V3 targets the RG extensor IN.',17)
    svg+=text(700,1470,'Other homonymous II and toe routes belong to motor/reflex circuits, not this central panel.',17)
    svg+='</svg>';save('W2L108_'+side+'_central_feedback',svg)
    table(DEST/('W2L108_'+side+'_central_feedback_edges.csv'),rows+feed)
    # A separate legible named-type key keeps conductance labels off the wires.
    key=start(1260,990)+text(630,42,f'W2L108 {side} · central feedback connection key',25)
    key+=text(30,80,'g: maximum conductance (uS); Erev, threshold, saturation: mV. Full IDs in CSV/SVG.',16,'start')
    headers=[('Key',30),('Source → target',105),('Named type',600),('g',930),('Erev',1010),('lo',1090),('hi',1170)]
    for title,x in headers:key+=text(x,122,title,16,'start')
    for i,e in enumerate(rows):
        y=162+i*43
        if i%2==0:key+=f'<rect x="20" y="{y-24}" width="1215" height="42" fill="#F1F4F5"/>'
        key+=text(30,y,e['figure_key'],16,'start')+text(105,y,e['source']+' → '+e['target'],15,'start')+text(600,y,e['type_name'],15,'start')
        for field,x in [('gmax_uS',930),('Erev_mV',1010),('threshold_mV',1090),('saturation_mV',1170)]:key+=text(x,y,f'{e[field]:.4g}',15,'start')
    key+=text(30,935,'Source: W2L108 verified standalone; receptor conversion and MN wiring are separate layers.',16,'start')+'</svg>'
    save('W2L108_'+side+'_central_feedback_key',key)

index={'source':str(SOURCE),'source_sha256':digest(SOURCE),'primary_figures':['W2L108_L_central_feedback','W2L108_R_central_feedback'],'source_neurons':108,'source_chemical_edges':216,'drawing_scope':'Direct Ia/II/Ib/contact projections onto six central half-centres per side, plus source RG-to-PF edges; internal circuitry grouped.','circuits_executed':0}
(DEST/'CENTRAL_FEEDBACK_DRAWING_RECORD.json').write_text(json.dumps(index,indent=2))
print(json.dumps(index))
