import io, sys
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

path = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\walk tester rearranged.aproj'
tree = ET.parse(path)
root = tree.getroot()

def walk(node, depth=0, maxdepth=4):
    cls = node.get('Type', '')
    nm_el = node.find('Name')
    nm = nm_el.text if nm_el is not None else node.get('Name', '')
    print('  ' * depth + f'<{node.tag}> "{nm}" {cls}')
    if depth < maxdepth:
        for ch in node:
            walk(ch, depth + 1, maxdepth)

for el in root.iter():
    if el.get('Type') == 'AnimatSim.Environment.NeuralModule':
        walk(el)
        break
