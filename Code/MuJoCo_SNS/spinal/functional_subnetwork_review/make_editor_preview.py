"""Build an isolated, source-preserving schematic-editor preview from saved data."""
import copy,json,hashlib
from pathlib import Path
import draw_full_motor_loops as F
import draw_feedback_v3 as V
from audit_saved_graphs import ROOT,digest

OUT=Path('D:/Users/Ben Bolen/NeuroGraphics_20261102/editor_preview')
OUT.mkdir(parents=True,exist_ok=True)
def rawdoc(title,source,raw,view):
    return {'schema':'connectome-schematic-views/1','identity':{'title':title,'source_file':str(source),'source_sha256':digest(source)},'sourceGraph':raw,'views':[view]}
def port(n,p):
    x,y,w,h=[n[k] for k in ('x','y','w','h')]
    choices=[('left',abs(p[0]-x),(p[1]-y)/h),('right',abs(p[0]-x-w),(p[1]-y)/h),('top',abs(p[1]-y),(p[0]-x)/w),('bottom',abs(p[1]-y-h),(p[0]-x)/w)]
    side,_,offset=min(choices,key=lambda z:z[1]);return {'side':side,'offset':offset}
def fixture(kind):
    p,es,d=V.load('w2lvar');ix=F.M.definitions(d)
    old=F.FullDrawing.output;F.FullDrawing.output=lambda self,*args:None
    try:slug,c=F.figure(ix,kind)
    finally:F.FullDrawing.output=old
    glyphs=[]
    for name,n in c.nodes.items():
        nn=copy.deepcopy(n);nn['id']=nn.pop('name');nn['nodeKeys']=c.aliases.get(name,[name] if any(x['name']==name for x in d['nodes']) else [])
        nn['shape']=c.icon_nodes.get(name,nn.get('shape','box'));glyphs.append(nn)
    wires=[]
    for i,w in enumerate(c.wires):
        def glyph_for(p,preferred):
            match=next((n for n in glyphs if preferred in n['nodeKeys'] and F.M.on_border(p,n)),None)
            return match or next(n for n in glyphs if F.M.on_border(p,n))
        a=glyph_for(w['points'][0],w['edges'][0]['source']);b=glyph_for(w['points'][-1],w['edges'][0]['target'])
        wires.append({'id':'wire:'+str(i),'from':a['id'],'to':b['id'],'fromPort':port(a,w['points'][0]),'toPort':port(b,w['points'][-1]),'bends':[{'x':x,'y':y} for x,y in w['points'][1:-1]],'records':copy.deepcopy(w['edges']),'sign':w['sign'],'color':w['color'],'dashed':w['dashed']})
    view={'title':c.title,'width':F.WIDTH,'height':c.h,'glyphs':glyphs,'wires':wires,'frames':[{'label':label,'x':x,'y':y,'w':w,'h':h} for label,x,y,w,h in c.groups],'annotations':[{'x':x,'y':y,'text':s,'size':sz} for x,y,s,sz in c.labels],'scope':'Historical saved w2lvar standing reference. Defined PF→MN routes have effective g=0. Layout is separate from source graph.'}
    return rawdoc('Saved w2lvar · '+('biarticular full loops' if kind=='biarticular' else 'shared PF and full muscle maps'),p,d,view)
def node(name,label,x,y,w=180,h=75,color='#9AD8F1',shape='box'):
    return {'id':name,'nodeKeys':[name],'label':label,'x':x,'y':y,'w':w,'h':h,'color':color,'shape':shape}
def wire(raw,i,a,b,pts,sign):
    return {'id':'wire:'+str(i),'from':a['id'],'to':b['id'],'fromPort':port(a,pts[0]),'toPort':port(b,pts[-1]),'bends':[{'x':x,'y':y} for x,y in pts[1:-1]],'records':[raw],'sign':sign,'color':'#333','dashed':False}
def current_rf():
    p=OUT/'sources/RFHAM_EXECUTED_SUBNETWORK_LEDGER.json';d=json.loads(p.read_text());ns=[]
    for stem,x in [('HIP',160),('KNEE',1030)]:
        for k,dx in [('E',0),('F',350)]:ns.append(node(f'PF_{stem}-{k}_l','ext' if k=='E' else 'flx',x+dx,160,180,80,'#FFE52E'))
    ns += [node('PF_BI_JOINT_RECT_l','RF PF',440,530,180,85,'#FFE52E'),node('PF_BI_JOINT_HAM_l','HAM PF',1110,530,180,85,'#FFE52E'),node('PF_IN_BI_PAIR_RECT_l','RF IN',760,650),node('PF_IN_BI_PAIR_HAM_l','HAM IN',760,800)]
    for name,label,x,y in [('rect_fem','Rectus femoris',250,1080),('semimem','Semimembranosus',1037.5,1030),('semiten','Semitendinosus',1250,1170),('bifemlh','Biceps femoris LH',1510,1030)]:ns.append(node('MN_'+name+'_l',label+'\nMN',x,y,225,85))
    by={n['id']:n for n in ns};ws=[]
    for i,e in enumerate(x for x in d['edges'] if x['side']=='l'):
        a,b=by[e['source']],by[e['target']];ax=a['x']+a['w']/2;bx=b['x']+b['w']/2
        if e['source'].startswith(('PF_HIP','PF_KNEE')):
            yy=330+i*34;pts=[(ax,a['y']+a['h']),(ax,yy),(bx,yy),(bx,b['y'])]
        elif e['target'].startswith('PF_IN'):
            yy=a['y']+a['h']/2;targety=b['y']+b['h']/2
            pts=[(a['x']+a['w'] if ax<bx else a['x'],yy),(bx,yy),(bx,b['y'])] if ax<bx else [(a['x'],yy),(1020,yy),(1020,targety),(b['x']+b['w'],targety)]
        elif e['source'].startswith('PF_IN'):
            yy=a['y']+a['h']/2;ty=b['y']+b['h']/2
            lane=1060 if bx>ax else 680;pts=[(a['x']+a['w'] if bx>ax else a['x'],yy),(lane,yy),(lane,ty),(b['x'] if bx>ax else b['x']+b['w'],ty)]
        else:
            if 'RECT' in e['source']:yy=960;pts=[(ax,a['y']+a['h']),(ax,yy),(bx,yy),(bx,b['y'])]
            else:
                sx={'MN_semimem_l':1150,'MN_semiten_l':1200,'MN_bifemlh_l':1250}[e['target']];yy={'MN_semimem_l':900,'MN_semiten_l':960,'MN_bifemlh_l':930}[e['target']]
                pts=[(sx,a['y']+a['h']),(sx,yy),(bx,yy),(bx,b['y'])]
        ws.append(wire(e,i,a,b,pts,'exc' if e['reversal_mV']>0 else 'inh'))
    v={'title':'Current selected human RF/HAM four-cell module','width':1850,'height':1370,'glyphs':ns,'wires':ws,'frames':[{'label':'Hip PF','x':100,'y':80,'w':650,'h':200},{'label':'K&A PF','x':970,'y':80,'w':650,'h':200}],'annotations':[],'scope':'Exact current selected left RF/HAM subnetwork: 12 cells, 12 synapses. Right homologues and full precision source records retained in sourceGraph. Sartorius/gracilis are separate.'}
    return rawdoc('Current selected RF/HAM · left module',p,d,v)
def current_syn6():
    p=OUT/'sources/Syn6_FINAL_COMPILED_GRAPH.json';d=json.loads(p.read_text());assert digest(p)=='697678b3281d7fb3b2a9322c2bc77b0671e7d05beb932fcb4364b8aa5ddaddea'
    ns=[node('HEEL_r','Heel contact',120,370,240,85,'#FFF0B3'),node('RG_E_r','ext HC',650,180,180,85,'#FFE52E'),node('RG_F_r','flx HC',1140,180,180,85,'#FFE52E'),node('InE_r','ext IN',650,370),node('InF_r','flx IN',1140,370)]
    by={n['id']:n for n in ns};ws=[]
    for a,b,pts in [('HEEL_r','InE_r',[(360,412.5),(650,412.5)]),('InE_r','RG_F_r',[(830,407.5),(1010,407.5),(1010,222.5),(1140,222.5)])]:
        e=next(e for e in d['chemical_edges'] if e['source']==a and e['target']==b);ws.append(wire(e,len(ws),by[a],by[b],pts,'exc' if e['reversal_mV']>0 else 'inh'))
    assert not any(e['source']=='HEEL_r' and e['target']=='InF_r' for e in d['chemical_edges'])
    v={'title':'Reviewed Syn6 · heel into extensor RG interneuron','width':1600,'height':740,'glyphs':ns,'wires':ws,'frames':[{'label':'RG','x':600,'y':100,'w':770,'h':410}],'annotations':[{'x':800,'y':600,'text':'Heel excites ipsilateral InE; InE inhibits RG-F. Internal RG motif connections are omitted from this focused view.','size':19}],'scope':'Current reviewed Syn6 composition: 794 cells, 3,200 active pairs. This view selects two pathways; full graph retained. No body/training-success claim.'}
    return rawdoc('Reviewed Syn6 · heel → InE → RG-F',p,d,v)

def main():
    docs=[fixture('multiple'),fixture('biarticular'),current_rf(),current_syn6()]
    manifest=[]
    for i,d in enumerate(docs):
        path=OUT/f'example_{i+1}.json';path.write_text(json.dumps(d,ensure_ascii=False),encoding='utf-8')
        manifest.append({'title':d['identity']['title'],'file':path.name,'source_sha256':d['identity']['source_sha256'],'view_glyphs':len(d['views'][0]['glyphs']),'view_wires':len(d['views'][0]['wires'])})
    library=json.dumps(docs,ensure_ascii=False).replace('<','\\u003c')
    template=(Path(__file__).parent/'editor_preview_shell.html').read_text(encoding='utf-8')
    (OUT/'Connectome_Schematic_Preview.html').write_text(template.replace('/*LIBRARY_JSON*/',library),encoding='utf-8')
    (OUT/'schematic_editor.js').write_bytes((Path(__file__).parent/'schematic_editor.js').read_bytes())
    (OUT/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(manifest))
if __name__=='__main__':main()
