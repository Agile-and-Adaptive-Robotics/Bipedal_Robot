import unittest,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from controller_library import standing,native108,reviewed_syn6,hierarchy
from project import flatten,normalize

class LibraryTests(unittest.TestCase):
    def test_all_four_full_graphs_survive_hierarchy(self):
        for family,flat,count in [('s3k',standing('s3k'),(410,1186)),('w2lvar',standing('w2lvar'),(888,7866)),('W2L108',native108(),(108,216)),('syn6',reviewed_syn6(),(794,3200))]:
            with self.subTest(family=family):
                saved=normalize(json.loads(json.dumps(hierarchy(flat))));nodes,edges,_=flatten(saved)
                self.assertEqual((len(nodes),len(edges)),count)
                self.assertEqual({n['id']:n['source_record'] for n in nodes},{n['id']:n['source_record'] for n in flat['nodes']})
                self.assertEqual({e['id']:e['source_record'] for e in edges},{e['id']:e['source_record'] for e in flat['edges']})
                self.assertEqual(saved['source'],flat['source'])
                self.assertEqual(sum(e['synapse']['gmax_uS']==0 for e in edges),sum(e['synapse']['gmax_uS']==0 for e in flat['edges']))

if __name__=='__main__':unittest.main()
