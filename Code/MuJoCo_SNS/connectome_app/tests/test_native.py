"""Static native-file checks only. No simulator is started."""
import copy,json,sys,tempfile,unittest,uuid
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import native_w2l as n
from project import FileSession,ProjectError,flatten
from graph_ops import flat_document

SOURCE=Path(__file__).resolve().parents[4]/'tmp/functional_subnetwork_remote/W2L108_verified/editable/W2L_VirtualWalker_Newer.aproj'
@unittest.skipUnless(SOURCE.exists(),'Native reference is machine-local')
class NativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.doc=n.read_model(SOURCE)
    def test_unchanged_bytes_and_actual_parameters(self):
        d=copy.deepcopy(self.doc);a,p=n.build_bytes(d)
        self.assertEqual(a,n.dec(d['native']['asim_b64']));self.assertEqual(p,SOURCE.read_bytes())
        ns,es,_=flatten(d);self.assertEqual((len(ns),len(es)),(108,216))
        self.assertEqual(next(x for x in ns if x['label']=='L Hip MN ext')['model']['rest_mV'],-100)
    def test_one_synapse_shared_type_clone_and_scalars(self):
        d=flat_document(self.doc);before=n.parse(n.dec(d['native']['asim_b64']));e=d['edges'][0];ident=e['id'];values=dict(gmax_uS=.12345678901234566,reversal_mV=-21.25,threshold_mV=-65.5,saturation_mV=-33.25);e['synapse']=values
        cell=d['nodes'][0];cell['model']['native_parameters']['tau_ms']=12.345;d['native']['adapters'][0]['gain']['C']=1.2345e-10
        a,p=n.build_bytes(d);r=n.parse(a);pr=n.parse(p);old={x.findtext('ID'):x for x in before.iter('Connexion')};now={x.findtext('ID'):x for x in r.iter('Connexion')};types={x.findtext('ID'):x for x in r.iter('SynapseType')}
        self.assertNotEqual(now[ident].findtext('SynapseTypeID'),old[ident].findtext('SynapseTypeID'))
        for k in old:
            if k!=ident:self.assertEqual(n.E.tostring(old[k]),n.E.tostring(now[k]))
        t=types[now[ident].findtext('SynapseTypeID')]
        self.assertEqual({k:n.number(t,v) for k,v in n.SYN.items()},values)
        gt=next(x for x in pr.iter('Link') if x.findtext('ID')==t.findtext('ID'))
        self.assertEqual(float(gt.find('MaxSynapticConductance').get('Actual')),values['gmax_uS']*1e-6)
        self.assertEqual(n.number(n.catalog(r)[cell['id']],'TimeConst'),12.345)
        # Entire rigid-body mechanics is preserved.
        self.assertEqual([n.E.tostring(x) for x in before.iter('RigidBody')],[n.E.tostring(x) for x in r.iter('RigidBody')])
    def test_add_rewire_delete_and_addflow_indices(self):
        d=flat_document(self.doc);first=copy.deepcopy(d['nodes'][0]);first['id']=str(uuid.uuid4());first['label']='New editable neuron';d['nodes'].append(first)
        edge=copy.deepcopy(d['edges'][0]);edge['id']=str(uuid.uuid4());edge['fromId']=first['id'];edge['from']=first['label'];d['edges'].append(edge)
        d['edges'][1]['toId']=first['id'];d['edges'][1]['to']=first['label']
        a,p=n.build_bytes(d);r=n.parse(a);pr=n.parse(p);self.assertEqual(len(n.catalog(r)),109);self.assertEqual(len(list(r.iter('Connexion'))),217)
        real=next(x for x in pr.iter('Node') if x.findtext('ID')==first['id'])
        for gate in ('CaActivation','CaDeactivation'):self.assertEqual(real.findtext(gate+'/ID'),n.catalog(r)[first['id']].findtext(gate+'/ID'))
        for owner in pr.iter():
            payload=owner.find('DiagramXml')
            if payload is None or not payload.text:continue
            flow=n.E.fromstring(payload.text.encode()).find('.//AddFlow')
            if flow is None:continue
            objects=[x for x in flow if x.tag in ('Node','Link')];by={x.findtext('ID'):x for x in owner.findall('Links/Link')}
            for item in objects:
                if item.tag!='Link':continue
                link=by[item.findtext('Tag')]
                for attr,endpoint in [('Org','OriginID'),('Dst','DestinationID')]:self.assertEqual(objects[int(item.get(attr))].findtext('Tag'),link.findtext(endpoint))
        d['nodes']=[x for x in d['nodes'] if x['id']!=first['id']];d['edges']=[e for e in d['edges'] if first['id'] not in(e['fromId'],e['toId'])];n.build_bytes(d)
    def test_save_as_companions_reopen_actual_truth_and_conflict(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'copy.aproj';s=FileSession();_,doc=s.save(self.doc,'a',p)
            self.assertTrue(p.with_suffix('.asim').exists());self.assertTrue(p.with_suffix('.editor.json').exists())
            self.assertEqual(SOURCE.read_bytes(),n.dec(self.doc['native']['project_b64']))
            reopened=s.open(p,'a');self.assertEqual(len(flatten(reopened)[0]),108)
            flat=flat_document(reopened);edgeid=flat['edges'][0]['id'];flat['edges'][0]['synapse']['gmax_uS']=.4321;s.save(flat,'a');again=s.open(p,'a');self.assertEqual(next(e for e in flatten(again)[1] if e['id']==edgeid)['synapse']['gmax_uS'],.4321)
            old=p.read_bytes();a=p.with_suffix('.asim');a.write_bytes(a.read_bytes()+b'\n')
            with self.assertRaises(ProjectError):s.save(again,'a')
            self.assertEqual(p.read_bytes(),old)
    def test_layout_and_ports_survive_native_save(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'layout.aproj';d=flat_document(self.doc);e=d['edges'][0];e['pts']=[{'x':12.5,'y':100}];e['fromPort']={'side':'bottom','offset':.234};e['toPort']={'angle':1.23}
            s=FileSession();s.save(d,'a',p);again=flat_document(s.open(p,'a'));actual=next(x for x in again['edges'] if x['id']==e['id'])
            for k in ('pts','fromPort','toPort'):self.assertEqual(actual[k],e[k])
    def test_removed_adapter_updates_both_formats(self):
        d=copy.deepcopy(self.doc);removed=d['native']['adapters'].pop(0)['id'];a,p=n.build_bytes(d)
        self.assertNotIn(removed,{x.findtext('ID') for x in n.parse(a).iter('Adapter')});pr=n.parse(p)
        self.assertNotIn(removed,{x.findtext('ID') for x in pr.iter('Node')})
        self.assertFalse(any(removed in(x.findtext('OriginID'),x.findtext('DestinationID')) for x in pr.iter('Link')))
if __name__=='__main__':unittest.main()
