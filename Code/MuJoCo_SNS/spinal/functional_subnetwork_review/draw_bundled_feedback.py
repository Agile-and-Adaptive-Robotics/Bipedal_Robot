"""Draw saved central sensory edges and explicit relays; never run a circuit."""
import json,html
from collections import defaultdict
from pathlib import Path
import pymupdf
from audit_saved_graphs import ROOT,OUT,table,digest
DEST=OUT/'drawings'
def txt(x,y,s,size=18,anchor='middle'):
    return ''.join(f'<text x="{x}" y="{y+i*size*1.25}" font-family="Arial" font-size="{size}" text-anchor="{anchor}">{html.escape(t)}</text>' for i,t in enumerate(s.split('\n')))
def box(x,y,w,h,color):return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{color}" stroke="#222" stroke-width="1.5"/>'
def save(name,svg):
    (DEST/(name+'.svg')).write_text(svg,encoding='utf-8');s=pymupdf.open(stream=svg.encode(),filetype='svg');d=pymupdf.open('pdf',s.convert_to_pdf());d.save(DEST/(name+'.pdf'));d[0].get_pixmap(matrix=pymupdf.Matrix(1.2,1.2)).save(DEST/(name+'.png'))
    d.close();s.close()
for family,scope in (('s3k','RG_PF'),('w2lvar','RG'),('w2lvar','PF'),('syn6','RG'),('syn6','PF')):
    src=ROOT/('Neuromechanical_Models/Mujoco_SNS_models/Gait2392_standing_families/'+family+'_graph_ledger.json');d=json.loads(src.read_text())
    central={n['name'] for n in d['nodes'] if n.get('layer') in ('RG','PF') or n['name'].startswith(('RG_','PF_','InE_','InF_'))}
    incoming=[e for e in d['chemical_edges'] if e['target'] in central and e['target'].endswith('_l') and (e['source'].startswith(('Ia_','II_','Ib_','HEEL_','TOEDF_','TOE_')))]
    if scope!='RG_PF':incoming=[e for e in incoming if e['target'].startswith('PF')==(scope=='PF')]
    groups=defaultdict(list)
    for e in incoming:
        channel='Ia' if e['source'].startswith('Ia_') else 'II' if e['source'].startswith('II_') else 'Ib' if e['source'].startswith('Ib_') else 'heel' if e['source'].startswith('HEEL_') else 'toe'
        groups[channel,e['target']].append(e)
    targets=sorted({e['target'] for e in incoming},key=lambda n:(0 if n.startswith('RG') else 1,n))
    extra=[]
    for target in list(targets):
        if 'IN' in target or target.startswith(('InE','InF')):
            extra.extend(e for e in d['chemical_edges'] if e['source']==target and e['target'] in central and e['target'].endswith('_l'))
    destinations=sorted({e['target'] for e in extra})
    height=max(1200,350+len(targets)*180)
    svg=f'<svg xmlns="http://www.w3.org/2000/svg" width="1500" height="{height}" viewBox="0 0 1500 {height}"><rect width="100%" height="100%" fill="white"/>'
    svg+=txt(750,38,f'{family} · left sensory feedback into {scope.replace("_"," / ")}',27)
    svg+=txt(750,73,'Muscle bundles are labelled by channel count; actual receiving interneurons stay explicit.',19)
    positions={k:(40,145+i*170) for i,k in enumerate(('Ia','II','Ib','heel','toe'))}
    target_pos={name:(770,130+i*180) for i,name in enumerate(targets)}
    tips=[];records=[];target_slots={}
    for i,((channel,target),es) in enumerate(groups.items()):
        x,y=positions[channel];tx,ty=target_pos[target]
        slot=target_slots.get(target,0);target_slots[target]=slot+1;ty+=18+slot*23
        lane=445+i*18
        # Toe has a real intervening TOEDF relay; do not bypass it.
        sx=x+275
        if channel=='toe':sx=430
        svg+=f'<polyline points="{sx},{y+40} {lane},{y+40} {lane},{ty} {tx},{ty}" fill="none" stroke="#333" stroke-width="{2.7 if len(es)>1 else 1.6}"/>'
        sign=es[0].get('sign_category','exc')
        tips.append(f'<polygon points="{tx-14},{ty} {tx},{ty-7} {tx},{ty+7}" fill="white" stroke="#222" stroke-width="1.5"/>' if sign=='exc' else f'<circle cx="{tx-6}" cy="{ty}" r="6" fill="black"/>')
        if len(es)>1:svg+=txt(tx-22,ty-6,f'{channel} ×{len(es)}',16,'end')
        record={'channel':channel,'target':target,'connection_count':len(es),'source_names':sorted({e['source'] for e in es}),'actual_edge_ids':[e['id'] for e in es],'full_edge_fields':es}
        records.append(record)
    for channel,(x,y) in positions.items():
        members={e['source'] for (ch,t),es in groups.items() if ch==channel for e in es}
        svg+=box(x,y,275,85,'#FFECA8' if members else '#EEEEEE')
        label=channel+'\n'+(f'{len(members)} muscle channels' if channel in ('Ia','II','Ib') and members else 'contact input' if members else 'No direct RG/PF route\nin this saved condition')
        svg+=txt(x+137.5,y+26,label,19)
        if channel=='toe' and members:
            relay='TOEDF_l';up=[e for e in d['chemical_edges'] if e['target']==relay]
            svg+=f'<line x1="315" y1="{y+40}" x2="330" y2="{y+40}" stroke="#333" stroke-width="1.6"/>'
            svg+=box(330,y+5,100,70,'#D3EABD')+txt(380,y+32,'TOEDF\nrelay / IN',16)
            svg+=f'<polygon points="316,{y+40} 330,{y+33} 330,{y+47}" fill="white" stroke="#222"/>'
            records.append({'channel':'toe relay upstream','full_edge_fields':up,'actual_edge_ids':[e['id'] for e in up]})
    for name,(x,y) in target_pos.items():
        inh='IN' in name or name.startswith(('InE','InF'))
        svg+=box(x,y,260,85,'#D3EABD' if inh else '#DDECF7' if name.startswith('RG') else '#FFE52E')+txt(x+130,y+33,name.replace('_l',' (L)').replace('_',' '),18)
        out=[e for e in extra if e['source']==name]
        for j,e in enumerate(out):
            dx=1180;dy=y+15+j*100
            svg+=f'<line x1="1030" y1="{y+40}" x2="{dx}" y2="{dy+30}" stroke="#333" stroke-width="1.6"/>'
            target_box=box(dx,dy,285,65,'#FFE52E')
            target_label=e['target'].replace('_l',' (L)').replace('_',' ')
            if e['target'] in targets:
                target_box=target_box.replace('stroke-width="1.5"','stroke-width="1.5" stroke-dasharray="6 4"')
                target_label+='\n(same HC; off-page reference)'
            svg+=target_box+txt(dx+142,dy+26,target_label,17)
            if e.get('sign_category')=='inh':svg+=f'<circle cx="{dx-6}" cy="{dy+30}" r="6" fill="black"/>'
            else:svg+=f'<polygon points="{dx-14},{dy+30} {dx},{dy+23} {dx},{dy+37}" fill="white" stroke="#222"/>'
            records.append({'channel':'interneuron output','full_edge_fields':[e],'actual_edge_ids':[e['id']]})
    svg+=''.join(tips)
    if not targets:svg+=box(760,300,580,190,'#F3F7EF')+txt(1050,355,'No direct Ia, II, Ib, heel, or toe\nchemical projection into RG/PF\nin this saved s3k standing condition.',23)
    footer=height-135
    if family=='syn6':
        svg+=txt(750,footer,'Saved heel/toe physical input gains are zero: the drawn chemical paths exist but receive no contact drive.',18)
    else:svg+=txt(750,footer,'A thick ×N line stands for N distinct source-muscle connections; the member list is in the companion key.',18)
    svg+=txt(750,footer+32,'Triangles have their base on the receiving object. Filled circles are inhibitory terminals.',18)
    svg+=txt(750,footer+64,'HC = half-centre · IN = interneuron · Dashed destination: reference to the same HC, not an added cell.',18)+'</svg>'
    save(family+'_central_feedback_'+scope,svg)
    (DEST/(family+'_central_feedback_'+scope+'_key.json')).write_text(json.dumps({'source':str(src),'sha256':digest(src),'bundles_and_relays':records},indent=2))
    if family=='w2lvar' and scope=='PF':
        roster=f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="1130"><rect width="100%" height="100%" fill="white"/>'+txt(600,45,'w2lvar · exact members of the muscle bundles',26)
        for column,channel in enumerate(('Ia','II','Ib')):
            members=sorted({e['source'] for (ch,target),es in groups.items() if ch==channel for e in es})
            x=40+column*390
            roster+=txt(x,100,f'{channel}: {len(members)} source channels',20,'start')
            for j,name in enumerate(members):roster+=txt(x,145+j*30,name,17,'start')
        roster+=txt(600,1060,'Each bundle expands into these separate muscle-channel connections, not one summed synapse.',18)
        roster+=txt(600,1090,'Endpoints, individual conductances, and edge IDs remain in the companion connection keys.',18)+'</svg>'
        save('w2lvar_muscle_bundle_members',roster)
print('Saved central-feedback bundle pages for s3k, w2lvar, and syn6; no circuits run.')
