"""Extract an editable working MN subnetwork without executing its source."""
import json,sys
from pathlib import Path
HERE=Path(__file__).parent;SPINAL=HERE.parent/'spinal'
sys.path.insert(0,str(SPINAL/'functional_subnetwork_review'))
import draw_full_motor_loops as F
import draw_feedback_v3 as V
from audit_saved_graphs import digest

def main():
    path,edges,data=V.load('w2lvar');ix=F.M.definitions(data)
    F.FullDrawing.output=lambda *args,**kwargs:None
    slug,c=F.figure(ix,'biarticular');raw_nodes={n['name']:n for n in data['nodes']}
    prefix_types={'Ia':'SN-Ia','II':'SN-II','Ib':'SN-Ib','IaIN':'IN-IaIN','IIX':'IN-IIe','IIIN':'IN-IIi','IBIN':'IN-IbIN','MN':'MN','RC':'RC'}
    nodes=[];ids={}
    for name,n in c.nodes.items():
        if name.startswith('PF_'):continue
        if name in raw_nodes:
            r=raw_nodes[name];p=r['source_neuron_parameters'];assert not p.get('gated') and not p.get('spiking')
            block=prefix_types[name.split('_')[0]];identifier=r['id'];model={'kind':'graded','capacitance_nF':p['membrane_capacitance']*1000,'conductance_uS':p['membrane_conductance'],'rest_mV':p['resting_potential'],'bias_nA':p['bias']}
            record={'source_record':r,'external_input':block.startswith('SN-') or block=='MN','external_output':block=='MN'}
        else:
            identifier='physical:'+name;block='MUSCLE' if name.startswith('muscle:') else 'SPINDLE' if name.startswith('spindle:') else 'TENDON' if name.startswith('tendon:') else 'ADAPTER';model={'kind':'physical'};record={}
        ids[name]=identifier
        nodes.append({'id':identifier,'type':block,'label':n.get('muscle_label') or n['label'] or ('Spindle' if block=='SPINDLE' else 'Force / Ib' if block=='TENDON' else block),'x':n['x']+n['w']/2,'y':n['y']+n['h']/2,'grp':'RF' if 'rect_fem' in name else 'SM','view':{'width':n['w'],'height':n['h'],'fill':n['color']},'model':model,**record})
    out=[]
    for w in c.wires:
        for e in w['edges']:
            if e['source'] not in ids or e['target'] not in ids:continue
            row={'id':e['id'],'fromId':ids[e['source']],'toId':ids[e['target']],'from':e['source'],'to':e['target'],'pts':[{'x':x,'y':y} for x,y in w['points'][1:-1]],'source_record':e,'tag':'','sign':w['sign'],'gain':e.get('gmax_uS',1)}
            # Explicit fractional boundary ports match the reference drawing.
            def port(name,p):
                n=c.nodes[name];x,y,ww,h=n['x'],n['y'],n['w'],n['h'];a=[('left',abs(p[0]-x),(p[1]-y)/h),('right',abs(p[0]-x-ww),(p[1]-y)/h),('top',abs(p[1]-y),(p[0]-x)/ww),('bottom',abs(p[1]-y-h),(p[0]-x)/ww)];side,_,off=min(a,key=lambda a:a[1]);return{'side':side,'offset':off}
            row['fromPort']=port(e['source'],w['points'][0]);row['toPort']=port(e['target'],w['points'][-1])
            if e.get('kind')=='physical_signal':row.update(kind='physical',map={'policy':e['policy'],'units':e['units']})
            else:row['synapse']={k:e[k] for k in ('gmax_uS','reversal_mV','threshold_mV','saturation_mV')}
            out.append(row)
    project={'schema':'sns-connectome-project/1','title':'Paired MN — editable working copy','clock':{'dt_ms':2.0},'presentation':{'style':'schematic'},'nodes':nodes,'edges':out,'groups':{'RF':'Rectus femoris local circuit','SM':'Semimembranosus local circuit'},'source':{'file':str(path),'sha256':digest(path),'condition':'Historical w2lvar standing local MN extraction; working copy, not the later selected Hunt controller.','units':'Source C numeric on seconds clock converted to nF for SNS millisecond clock: C_nF=1000*C_numeric.'},'interface':{'inputs':'Ia/II/Ib and direct MN currents; omitted PF/RG/other-muscle boundary drives are supplied externally.','physical_maps':data['encoder_units_and_signs'],'actuation':data['actuation']}}
    destination=HERE/'examples/paired_mn.sns.json';destination.parent.mkdir(exist_ok=True);destination.write_text(json.dumps(project,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({'neurons':sum(n['model']['kind']=='graded' for n in nodes),'physical_symbols':sum(n['model']['kind']=='physical' for n in nodes),'neural_synapses':sum(e.get('kind')!='physical' for e in out),'physical_maps':sum(e.get('kind')=='physical' for e in out)}))
if __name__=='__main__':main()
