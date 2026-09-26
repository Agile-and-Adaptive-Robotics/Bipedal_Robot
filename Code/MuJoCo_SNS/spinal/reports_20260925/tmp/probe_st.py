import io, sys
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
path = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\walk tester rearranged.aproj'
tree = ET.parse(path)
root = tree.getroot()

# print one SpikingChemical synapse-type Link raw
n = 0
for lk in root.iter('Link'):
    cls = lk.findtext('ClassName') or ''
    if 'SpikingChemical' in cls:
        print(ET.tostring(lk, encoding='unicode')[:2500])
        n += 1
        if n >= 2:
            break
# where does <Nodes> live? print path of parents
def parent_map(root):
    return {c: p for p in root.iter() for c in p}
pm = parent_map(root)
for el in root.iter('Nodes'):
    p = pm.get(el)
    gp = pm.get(p)
    print("\n<Nodes> parent:", p.tag if p is not None else None,
          "grandparent:", gp.tag if gp is not None else None)
