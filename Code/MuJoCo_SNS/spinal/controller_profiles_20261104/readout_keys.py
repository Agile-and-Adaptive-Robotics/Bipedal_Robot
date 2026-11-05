"""Stable mirrored keys for existing PF→MN connections; no new topology."""
from muscle_map import classify
def edge_key(source,target):
 act=target[3:];side=act[-1];root=source[:-2]if source.endswith(('_l','_r'))else source
 relation=('ipsi'if source[-1]==side else'contra')if source.endswith(('_l','_r'))else'mid'
 return classify(act).base+'|'+relation+'|'+root
def multiplier(source,target,saved):
 best=saved.get('best',saved)
 return best.get('gains',{}).get(classify(target[3:]).base,1.)*best.get('edge_gains',{}).get(edge_key(source,target),1.)
