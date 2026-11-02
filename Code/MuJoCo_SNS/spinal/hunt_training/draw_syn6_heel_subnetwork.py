"""Editable scoped neural schematic from the actual repaired compiled graph."""
from pathlib import Path
import argparse,csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Circle,FancyArrowPatch,PathPatch
from matplotlib.path import Path as DrawingPath
import matplotlib.patheffects as effects


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--graph',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();raw=args.graph.read_bytes();graph=json.loads(raw);out=args.out;out.mkdir(parents=True,exist_ok=True)
    wanted=('HEEL_r','RG_E_r','RG_F_r','InE_r','InF_r')
    nodes={n['name']:n for n in graph['nodes']if n['name']in wanted}
    if set(nodes)!=set(wanted):raise ValueError('Required actual source nodes are absent')
    edges={(e['source'],e['target']):e for e in graph['chemical_edges']if e['source']in nodes and e['target']in nodes}
    specs=[('S1','RG_E_r','RG_F_r',[(251,108),(273,87),(319,87),(341,108)],'curve',(296,81)),
           ('S2','RG_F_r','RG_E_r',[(341,122),(316,150),(276,150),(251,122)],'curve',(296,153)),
           ('S3','InE_r','RG_F_r',[(229,207),(351,137)],'line',(260,201)),
           ('S4','InF_r','RG_E_r',[(363,207),(240,137)],'line',(333,201)),
           ('S5','RG_E_r','InE_r',[(215,134),(215,202)],'line',(196,169)),
           ('S6','RG_F_r','InF_r',[(377,134),(377,202)],'line',(398,169)),
           ('S7','HEEL_r','InE_r',[(92.5,223),(194,223)],'line',(131,212))]
    if set(edges)!={(s,t)for _,s,t,*_ in specs}:raise ValueError('Actual graph differs from this scoped routing; do not invent or omit internal edges')
    heel_central=[e for e in graph['chemical_edges']if e['source']=='HEEL_r'and e['target'].startswith(('RG_','InE_','InF_'))]
    if[e['target']for e in heel_central]!=['InE_r']:raise ValueError('Source contains another heel-to-RG/IN branch')
    adapter=[a for a in graph['input_adapters']if a['target']=='HEEL_r']
    if len(adapter)!=1 or adapter[0]['input']!='HEEL_c_r'or adapter[0]['weight']!=1.:raise ValueError('Unexpected actual heel input map')
    plt.rcParams.update({'font.family':'Arial','font.size':10.5,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(6.5,4.5));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,468),ylim=(324,0));ax.axis('off')
    ax.text(234,20,'Repaired Syn6: right heel / RG subnetwork',ha='center',va='center',fontsize=12,fontweight='bold')
    ax.text(234,40,'Scoped view; PF branches and other inputs are omitted.',ha='center',va='center',fontsize=10.5)
    group=FancyBboxPatch((150,60),292,200,boxstyle='round,pad=0,rounding_size=12',facecolor='#F6F8FA',edgecolor='#B1BAC4',linewidth=1.)
    group.set_gid('right_RG_subnetwork');ax.add_patch(group)
    ax.text(165,80,'Right RG',fontsize=10.5,color='#46515B')
    ledger=[]
    for sid,source,target,points,shape,label in specs:
        e=edges[(source,target)];exc=e['reversal_mV']>nodes[target]['rest_mV']
        path=DrawingPath(points,[DrawingPath.MOVETO]+([DrawingPath.CURVE4]*3 if shape=='curve'else[DrawingPath.LINETO]))
        if exc:artist=FancyArrowPatch(path=path,arrowstyle='-|>',mutation_scale=12,linewidth=1.3,color='#222222',zorder=2)
        else:
            artist=PathPatch(path,fill=False,linewidth=1.3,color='#222222',zorder=2)
            if sid=='S4':artist.set_path_effects([effects.Stroke(linewidth=4.5,foreground='#F6F8FA'),effects.Normal()])
            endpoint=Circle(points[-1],radius=3.,facecolor='#222222',edgecolor='#222222',zorder=5);endpoint.set_gid('inhibitory_terminal_'+sid);ax.add_patch(endpoint)
        artist.set_gid('synapse_'+sid+'_'+source+'_to_'+target);ax.add_patch(artist)
        ax.text(*label,sid,ha='center',va='center',fontsize=10.5,bbox={'facecolor':'#F6F8FA'if label[0]>150 else'white','edgecolor':'none','pad':.5},zorder=6)
        ledger.append(dict(id=sid,**e,terminal='triangle'if exc else'filled_circle'))
    # Actual external-current adapter; this is not a chemical synapse.
    port=FancyBboxPatch((18,126),84,38,boxstyle='round,pad=0,rounding_size=4',facecolor='white',edgecolor='#687681',linewidth=1.)
    port.set_gid('actual_input_HEEL_c_r');ax.add_patch(port)
    ax.text(60,139,'HEEL_c_r',ha='center',va='center',fontsize=10.5);ax.text(60,153,'current (nA)',ha='center',va='center',fontsize=10.5)
    incoming=FancyArrowPatch((60,164),(60,205),arrowstyle='->',mutation_scale=10,linewidth=1.1,color='#687681',zorder=2)
    incoming.set_gid('input_adapter_HEEL_c_r_to_HEEL_r');ax.add_patch(incoming);ax.text(78,187,'×1',fontsize=10.5,ha='center',va='center')
    shapes=[('HEEL_r',(60,223),'HEEL','#E8EDF1','box'),('RG_E_r',(215,115),'RG-E','#DCEBF5','box'),
            ('RG_F_r',(377,115),'RG-F','#F5DDE0','box'),('InE_r',(215,223),'InE','#FFF1BD','circle'),('InF_r',(377,223),'InF','#FFF1BD','circle')]
    for name,(x,y),label,color,kind in shapes:
        if kind=='circle':patch=Circle((x,y),21,facecolor=color,edgecolor='#333333',linewidth=1.2,zorder=4)
        else:
            width,height=(65,36)if name=='HEEL_r'else(72,38)
            patch=FancyBboxPatch((x-width/2,y-height/2),width,height,boxstyle='round,pad=0,rounding_size=10',facecolor=color,edgecolor='#333333',linewidth=1.2,zorder=4)
        patch.set_gid('node_'+name);ax.add_patch(patch);ax.text(x,y,label,ha='center',va='center',fontsize=12,fontweight='bold',zorder=5)
    ax.add_patch(FancyArrowPatch((25,283),(55,283),arrowstyle='-|>',mutation_scale=12,color='#222222',linewidth=1.3));ax.text(67,283,'Excitatory',va='center',fontsize=10.5)
    ax.plot([170,200],[283,283],color='#222222',linewidth=1.3);ax.add_patch(Circle((200,283),3,color='#222222'));ax.text(212,283,'Inhibitory',va='center',fontsize=10.5)
    ax.add_patch(FancyArrowPatch((322,283),(350,283),arrowstyle='->',mutation_scale=10,color='#687681',linewidth=1.1));ax.text(362,283,'Input',va='center',fontsize=10.5)
    ax.text(234,309,'All cells graded; RG half-centers retain persistent-Na dynamics.',ha='center',va='center',fontsize=10.5)
    stem=out/'SYN6_RIGHT_HEEL_RG'
    fig.savefig(stem.with_suffix('.svg'),metadata={'Description':'Actual compiled five-node, seven-synapse induced subgraph; source IDs and exact fields in PROVENANCE.json.'})
    fig.savefig(stem.with_suffix('.pdf'));fig.savefig(stem.with_suffix('.png'),dpi=200);plt.close(fig)
    fields=['id','source','target','source_index','target_index','gmax_uS','reversal_mV','threshold_mV','saturation_mV','terminal']
    with(out/'SYNAPSE_FIELDS.csv').open('w',newline='',encoding='utf-8')as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(ledger)
    tablefig=plt.figure(figsize=(6.5,3.5));ta=tablefig.add_axes([0,0,1,1]);ta.set(xlim=(0,468),ylim=(252,0));ta.axis('off')
    ta.text(234,18,'Actual compiled synapse fields',ha='center',va='center',fontsize=12,fontweight='bold')
    widths=[22,65,65,115,57,57,57];left=15;xs=[]
    for width in widths:xs.append(left+width/2);left+=width
    headers=['ID','Source','Target','gmax / µS','Erev / mV','θ / mV','Sat / mV']
    for x,header in zip(xs,headers):ta.text(x,44,header,ha='center',va='center',fontsize=10.5,fontweight='bold')
    for i,e in enumerate(ledger):
        y=67+i*23
        values=[e['id'],e['source'],e['target'],repr(e['gmax_uS']),repr(e['reversal_mV']),repr(e['threshold_mV']),repr(e['saturation_mV'])]
        if i%2==0:ta.axhspan(y-10,y+10,color='#F3F5F7',zorder=0)
        for x,value in zip(xs,values):ta.text(x,y,value,ha='center',va='center',fontsize=10.5)
    ta.text(234,238,'Fields are unrounded, round-trip values from the executed graph.',ha='center',va='center',fontsize=10.5)
    tablefig.savefig(out/'SYNAPSE_FIELDS.pdf');tablefig.savefig(out/'SYNAPSE_FIELDS.png',dpi=200);plt.close(tablefig)
    md=['# Exact right heel/RG synapse fields','','| ID | Source | Target | gmax (µS) | Erev (mV) | Threshold (mV) | Saturation (mV) |','| --- | --- | --- | ---: | ---: | ---: | ---: |']
    for e in ledger:md.append('| '+' | '.join(str(e[k])for k in('id','source','target','gmax_uS','reversal_mV','threshold_mV','saturation_mV'))+' |')
    (out/'SYNAPSE_FIELDS.md').write_text('\n'.join(md),encoding='utf-8')
    external=[e for e in graph['chemical_edges']if(e['source']in nodes)!=(e['target']in nodes)]
    provenance=dict(source_graph=str(args.graph.resolve()),source_graph_sha256=hashlib.sha256(raw).hexdigest(),
                    source_profile_name=graph['name'],source_edge_operator_sha256=graph['edge_operator_sha256'],
                    nodes=[nodes[name]for name in wanted],chemical_edges=ledger,input_adapter=adapter[0],
                    all_internal_connections_shown=True,internal_edge_count=len(ledger),
                    external_connections_omitted=external,removed_heel_InF_branches=graph['removed_branches'],
                    scope='Right five-cell induced heel/RG subgraph; full PF and other external branches are not drawn. No dangling off-page wires.',
                    svg_text_editable=True,size_inches=[6.5,4.5],minimum_font_pt=10.5,
                    runtime_or_body_changes=False,generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    caption=('Repaired Syn6 right heel/rhythm-generation subnetwork. The actual heel-current input drives the graded HEEL encoder; its RG-interneuron branch excites InE only. InE inhibits RG-F, whereas the recurrent RG-F→InF→RG-E pathway is retained. Weak reciprocal RG excitation is also shown. Triangles indicate excitatory reversal relative to rest, and filled circles indicate inhibitory reversal; actual current remains g·s(Vpre)·(Erev−Vpost). All seven chemical connections among the five displayed cells are drawn from the final compiled graph. The retained heel→PF_IN_E branch, full PF circuitry, descending inputs, and commissural connections are outside this scoped view. Exact fields and source identities accompany the figure.\n')
    (out/'CAPTION.md').write_text(caption,encoding='utf-8')
    print(json.dumps(dict(output=str(out),nodes=len(nodes),chemical_edges=len(ledger),input_adapters=1,source_graph_sha256=provenance['source_graph_sha256'])))


if __name__=='__main__':main()
