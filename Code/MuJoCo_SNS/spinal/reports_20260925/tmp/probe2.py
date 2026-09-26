import io, sys
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
path = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\walk tester rearranged.aproj'
root = ET.parse(path).getroot()

def txt(el, tag):
    x = el.find(tag)
    return x.text if x is not None else None

ns = root.find('.//NervousSystem')
print("NervousSystem children:", [c.tag for c in ns])
node = ns.find('Node')
print("module Node ClassName:", txt(node, 'ClassName'), " ID:", txt(node, 'ID'))
print("module Node children:", [c.tag for c in node])
nodes = node.find('Nodes')
print("n <Node> in Nodes:", len(list(nodes)))
from collections import Counter
cc = Counter((txt(nd, 'ClassName') or '?').split('.')[-1] for nd in nodes)
print(cc)
