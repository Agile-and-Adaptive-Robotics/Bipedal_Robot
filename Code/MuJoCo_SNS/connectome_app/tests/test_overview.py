import copy,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from graph_ops import bilateral_overview,flat_document
from project import normalize,flatten

def legacy_overview():
    nodes=[];edges=[]
    for layer in ('RG','PF','MN','Afferents'):
        sides=[]
        for side in ('Left','Right'):
            key=layer+':'+side;leaf={'id':key,'type':'MN' if layer=='MN' else 'SN-Ia' if layer=='Afferents' else 'HC-'+layer+'-E','label':key,'layer':layer,'x':10,'y':10,'model':{'kind':'graded','rest_mV':0}}
            sides.append({'id':'side:'+key,'type':'SUB','label':side+' · 1 cells','scope':{'layer':layer,'side':side},'x':10,'y':10,'sub':{'nodes':[leaf],'edges':[]}})
        nodes.append({'id':'layer:'+layer,'type':'SUB','label':layer+'\n2 cells','scope':{'layer':layer},'x':700,'y':240 if layer=='RG' else 420 if layer=='PF' else 760 if layer=='MN' else 590,'sub':{'nodes':sides,'edges':[]}})
    for side in ('Left','Right'):
        edges.append({'id':'feedback:'+side,'fromId':'layer:Afferents','fromPath':['side:Afferents:'+side,'Afferents:'+side],'toId':'layer:RG','toPath':['side:RG:'+side,'RG:'+side],'synapse':{'gmax_uS':.5,'reversal_mV':8,'threshold_mV':0,'saturation_mV':5},'pts':[]})
    edges.append({'id':'cross','fromId':'layer:RG','fromPath':['side:RG:Left','RG:Left'],'toId':'layer:RG','toPath':['side:RG:Right','RG:Right'],'synapse':{'gmax_uS':.75,'reversal_mV':-5,'threshold_mV':0,'saturation_mV':5},'pts':[]})
    return normalize({'nodes':nodes,'edges':edges,'inputs':[{'name':'left current','node_id':'layer:RG','node_path':['side:RG:Left','RG:Left']}],'presentation':{'style':'schematic'}})

class OverviewTests(unittest.TestCase):
    def test_existing_overview_promotes_sides_and_retains_afferents(self):
        old=legacy_overview();prior=copy.deepcopy(old);new=bilateral_overview(old)
        self.assertEqual(old,prior)
        self.assertEqual([(n['scope']['layer'],n['scope']['side']) for n in new['nodes'] if n['scope']['layer']=='RG'],[('RG','Left'),('RG','Right')])
        self.assertEqual(len(new['nodes']),8)
        self.assertEqual({n['id']:n['model'] for n in flatten(new)[0]},{n['id']:n['model'] for n in flatten(old)[0]})
        before=flat_document(old);after=flat_document(new)
        def semantics(g):return {e['id']:(e['fromId'],e['toId'],e['synapse']) for e in g['edges']}
        self.assertEqual(semantics(before),semantics(after))
        for side in ('Left','Right'):
            e=next(e for e in new['edges'] if e['id']=='feedback:'+side)
            self.assertEqual(e['toId'],'side:RG:'+side);self.assertEqual(e['toPath'],['RG:'+side]);self.assertEqual(e['fromId'],'side:Afferents:'+side)
        self.assertEqual(new['inputs'][0]['node_id'],'side:RG:Left');self.assertEqual(new['inputs'][0]['node_path'],['RG:Left'])
        self.assertEqual(before['inputs'],after['inputs'])
        by={(n['scope']['layer'],n['scope']['side']):n for n in new['nodes']}
        for side in ('Left','Right'):self.assertLess(by['RG',side]['y'],by['PF',side]['y']);self.assertLess(by['PF',side]['y'],by['MN',side]['y'])
    def test_user_regrouping_after_migration_remains_intentional(self):
        new=bilateral_overview(legacy_overview());new['nodes'][0]['label']='My edited group'
        self.assertEqual(bilateral_overview(new),new)
    def test_plain_models_are_not_reorganized(self):
        doc=normalize({'nodes':[{'id':'a','type':'MN','label':'My cell'}],'edges':[]})
        self.assertEqual(bilateral_overview(doc),doc)
if __name__=='__main__':unittest.main()
