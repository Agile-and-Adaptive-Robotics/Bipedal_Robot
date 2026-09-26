import io, sys, re
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

path = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\walk tester rearranged.aproj'
tree = ET.parse(path)
root = tree.getroot()

def actual(el):
    """Return Actual value of a Value/Scale/Actual triplet element."""
    if el is None:
        return None
    a = el.find('Actual')
    return a.text if a is not None else el.text

# ---- enumerate all Environment + Organism + NeuralModule structure at high level
def walk(node, depth=0, maxdepth=3):
    tag = node.tag
    cls = node.get('Type', '')
    nm = node.find('Name')
    nm = nm.text if nm is not None else node.get('Name', '')
    print('  ' * depth + f'<{tag}> {nm} {cls}')
    if depth < maxdepth:
        for ch in node:
            walk(ch, depth + 1, maxdepth)

# find the nervous system node
for el in root.iter():
    if el.tag == 'Node' and el.get('Type') == 'AnimatSim.Environment.NeuralModule':
        walk(el, 0, 2)
        break
