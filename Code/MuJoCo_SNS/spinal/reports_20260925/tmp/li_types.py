import io, sys
import xml.etree.ElementTree as ET
from collections import Counter
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

path = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\walk tester rearranged.aproj'
tree = ET.parse(path)
root = tree.getroot()
c = Counter(el.tag for el in root.iter())
print("tags:", dict(c))
tc = Counter(el.get('Type') for el in root.iter() if el.get('Type'))
for k, v in sorted(tc.items()):
    print(f"{v:5d}  {k}")
