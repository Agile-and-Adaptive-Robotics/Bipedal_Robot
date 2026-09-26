import io, sys
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

path = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\walk tester rearranged.aproj'
tree = ET.parse(path)
root = tree.getroot()

def txt(el, tag):
    x = el.find(tag)
    return x.text if x is not None else None

def actual(el, tag):
    x = el.find(tag)
    if x is None:
        return None
    a = x.find('Actual')
    return a.text if a is not None else x.text

# ---- NervousSystem module-level params
ns = None
for el in root.iter('NervousSystem'):
    ns = el
    break
print("=== NervousSystem ===")
for ch in ns:
    if ch.tag in ('TimeStep', 'AHPEquilibriumPotential', 'SpikePeak', 'SpikeStrength',
                  'CaEquilibriumPotential', 'RefractoryPeriod', 'UseCriticalPeriod',
                  'StartCriticalPeriod', 'EndCriticalPeriod', 'TTX', 'Cd', 'HH'):
        print(f"  {ch.tag} = {ch.text}")

# ---- neurons (Node elements inside Diagram Nodes)
print("\n=== NEURONS (Node) ===")
for el in root.iter('Node'):
    nm = txt(el, 'NodeName') or ''
    # Node elements: find class
    cls = None
    for c in el.iter('ClassName'):
        cls = c.text
        break
    if cls is None:
        # try PartType
        pt = el.find('.//PartType')
        cls = pt.text if pt is not None else '?'
    # gather all params
    vals = {}
    for tag in ('RestingPotential', 'TimeConstant', 'InitialThreshold', 'RelativeAccomodation',
                'AccomodationTimeConstant', 'AHP_Conductance', 'AHP_TimeConstant',
                'MaxCaConductance', 'TonicStimulus', 'RelativeSize', 'Cd', 'HH'):
        v = actual(el, tag)
        if v is not None:
            vals[tag] = v
    name_el = None
    # neuron display name: <Node> block contains Text? or a Name child
    for t in el.iter('Text'):
        name_el = t.text
        break
    print(f"  [{el.get('ID')}] class={cls} text={name_el!r} params={vals}")

print("\n=== SYNAPSE TYPES ===")
stypes = {}
for el in root.iter():
    if el.tag == 'SynapseType' or (el.tag == 'Simulation' and el.get('Type') and 'SynapseType' in str(el.get('Type'))):
        pass
# synapse types may be <Simulation Type="...SynapseType"> blocks
for el in root.iter('Simulation'):
    t = el.get('Type') or ''
    if 'SynapseType' in t:
        eq = txt(el, 'EquilibriumPotential')
        print(f"  {t}  EquilibriumPotential={eq}")

print("\n=== LINKS (connexions) ===")
for el in root.iter('Link'):
    sid = txt(el, 'OriginID')
    did = txt(el, 'DestinationID')
    st = txt(el, 'SynapticTypeID')
    g = actual(el, 'SynapticConductance')
    d = txt(el, 'ConductionDelay')
    print(f"  {sid} -> {did}  type={st} g={g} delay={d}")

print("\n=== STIMULUS ===")
for el in root.iter('Stimulus'):
    for tag in ('Name', 'StartTime', 'EndTime', 'AlwaysActive', 'ValueType', 'Equation',
                'Current', 'PositionX', 'PositionY', 'PositionZ'):
        v = txt(el, tag)
        if v is not None:
            print(f"  {tag} = {v}")

print("\n=== ADAPTER LINKS (DataTypeID...) ===")
for el in root.iter('Adapter'):
    print(" adapter:", ET.tostring(el, encoding='unicode')[:600])

print("\n=== RECEPTIVE FIELD SENSOR ===")
for el in root.iter('ReceptiveFieldSensor'):
    print(" ", ET.tostring(el, encoding='unicode')[:800])

# neuron IDs -> names: the <Node ID=..> has child <Name>? let's print raw small block for first node
print("\n=== first Node raw ===")
for el in root.iter('Node'):
    s = ET.tostring(el, encoding='unicode')
    print(s[:1500])
    break
