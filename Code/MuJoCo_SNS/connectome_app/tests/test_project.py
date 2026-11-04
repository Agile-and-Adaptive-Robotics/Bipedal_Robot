import copy,json,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from project import FileSession,ProjectError,normalize,flatten
from runtime import build_project
FIXTURE=Path(__file__).resolve().parents[1]/'examples/paired_mn.sns.json'

class ProjectTests(unittest.TestCase):
    def setUp(self):self.doc=normalize(json.loads(FIXTURE.read_text(encoding='utf-8')))
    def test_roundtrip_complete_parameters_and_copy(self):
        with tempfile.TemporaryDirectory() as folder:
            a=Path(folder)/'original.json';b=Path(folder)/'copy.json';session=FileSession();session.save(self.doc,'a',a);original=a.read_bytes();opened=session.open(a,'a')
            edge=next(e for e in opened['edges'] if e.get('kind')!='physical');edge['synapse']['gmax_uS']=0.12345678901234566
            saved,doc=session.save(opened,'a',b);self.assertEqual(a.read_bytes(),original);self.assertEqual(session.open(b,'b'),doc)
            self.assertEqual(next(e for e in doc['edges'] if e['id']==edge['id'])['source_record'],edge['source_record'])
    def test_external_edit_not_overwritten(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'a.json';s=FileSession();s.save(self.doc,'a',p);p.write_text('{}')
            with self.assertRaises(ProjectError):s.save(self.doc,'a')
            self.assertEqual(p.read_text(),'{}')
    def test_dangling_and_compiled_import_rejected(self):
        bad=copy.deepcopy(self.doc);bad['edges'][0]['toId']='missing'
        with self.assertRaises(ProjectError):normalize(bad)
        with self.assertRaises(ProjectError):normalize({'nodes':[],'chemical_edges':[]})
    def test_saved_graph_is_consumed_by_builder(self):
        edge=next(e for e in self.doc['edges'] if e.get('kind')!='physical');edge['synapse']['gmax_uS']=0.12345678901234566
        model=build_project(self.doc)
        self.assertEqual((model['neurons'],model['synapses'],len(model['output_names'])),(18,30,2))
        built=[e['params']['max_conductance'] for e in model['network'].connections]
        self.assertIn(edge['synapse']['gmax_uS'],built);self.assertEqual(model['steps_executed'],0)
    def test_nested_child_endpoint_is_buildable(self):
        doc=copy.deepcopy(self.doc);selected=doc['nodes'][:3];ids={n['id'] for n in selected};inner=[e for e in doc['edges'] if e['fromId'] in ids and e['toId'] in ids]
        host={'id':'module','type':'SUB','label':'Module','x':0,'y':0,'sub':{'nodes':selected,'edges':inner}}
        doc['nodes']=[n for n in doc['nodes'] if n['id'] not in ids]+[host];doc['edges']=[e for e in doc['edges'] if e not in inner]
        for e in doc['edges']:
            for side in ('from','to'):
                if e[side+'Id'] in ids:e[side+'Path']=[e[side+'Id']];e[side+'Id']='module'
        doc['inputs']=[{'name':'nested afferent','node_id':'module','node_path':[selected[0]['id']]}]
        model=build_project(doc);self.assertEqual((model['neurons'],model['synapses']),(18,30));self.assertEqual(model['input_names'],['nested afferent'])
    def test_unimplemented_model_is_not_silently_changed(self):
        n=next(n for n in self.doc['nodes'] if n['model']['kind']=='graded');n['model']['kind']='persistent_NaP'
        self.assertEqual(normalize(self.doc)['nodes'][0]['model']['kind'],'persistent_NaP')
        with self.assertRaises(ProjectError):build_project(self.doc)
if __name__=='__main__':unittest.main()
