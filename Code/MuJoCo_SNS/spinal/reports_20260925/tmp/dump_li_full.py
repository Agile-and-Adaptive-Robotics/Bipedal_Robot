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

# ---- find the IntegrateFire NeuralModule
nm = None
for el in root.iter('Node'):
    c = txt(el, 'ClassName')
    if c == 'IntegrateFireGUI.DataObjects.Behavior.NeuralModule':
        nm = el
        break
print("neural module timestep:", actual(nm, 'TimeStep'))

# ---- synapse types: Link elements under <SynapseTypes>
stypes = {}
for st in nm.find('SynapseTypes'):
    sid = txt(st, 'ID')
    cls = txt(st, 'ClassName')
    stypes[sid] = dict(
        cls=cls.split('.')[-1],
        eq=actual(st, 'EquilibriumPotential'),
        maxg=actual(st, 'MaxSynapticConductance'),
        thresh=actual(st, 'PreSynapticThreshold'),
        sat=actual(st, 'PreSynapticSaturationLevel'),
    )
print("\n=== SYNAPSE TYPES ===")
for k, v in stypes.items():
    print(f"  {k}  {v}")

# ---- neurons: Node elements with ClassName Spiking/NonSpiking under <Nodes>
print("\n=== NEURONS ===")
neurons = {}
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if 'Neurons.' not in cls:
        continue
    nid = txt(nd, 'ID')
    name = txt(nd, 'Text') or txt(nd, 'Name') or ''
    pars = dict(
        cls=cls.split('.')[-1],
        rest=actual(nd, 'RestingPotential'),
        tc=actual(nd, 'TimeConstant'),
        thr=actual(nd, 'InitialThreshold'),
        acc=actual(nd, 'RelativeAccomodation'),
        acctc=actual(nd, 'AccomodationTimeConstant'),
        ahp_g=actual(nd, 'AHP_Conductance'),
        ahp_tc=actual(nd, 'AHP_TimeConstant'),
        ca_g=actual(nd, 'MaxCaConductance'),
        tonic=actual(nd, 'TonicStimulus'),
    )
    # Ca activation/deactivation gates: MidPoint/Slope under CaActivation/CaDeactivation
    for gate in ('CaActivation', 'CaDeactivation'):
        g = nd.find(gate)
        if g is not None:
            pars[gate] = (actual(g, 'MidPoint'), actual(g, 'Slope'))
    neurons[nid] = dict(name=name, **pars)
    print(f"  {name!r:32s} {pars}")

# ---- connexions: Link elements under <Links>
print("\n=== CONNEXIONS ===")
nbyid = {k: v['name'] for k, v in neurons.items()}
for lk in nm.find('Links'):
    o = txt(lk, 'OriginID')
    d = txt(lk, 'DestinationID')
    st = txt(lk, 'SynapticTypeID')
    g = actual(lk, 'SynapticConductance')
    on, dn = nbyid.get(o, o), nbyid.get(d, d)
    stn = stypes.get(st, {}).get('cls', st)
    print(f"  {on!r} -> {dn!r}  type={stn} eq={stypes.get(st,{}).get('eq')} g={g}")

# ---- adapters: outside the neural module (organism level) — find all NodeToPhysicalAdapter / PhysicalToNodeAdapter
print("\n=== ADAPTERS ===")
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if 'Adapter' in cls:
        name = txt(nd, 'Text') or txt(nd, 'Name') or ''
        # subnodes hold the actual adapter config
        subs = []
        for sub in nd.find('Nodes') or []:
            scls = txt(sub, 'ClassName') or ''
            dt = txt(sub, 'DataTypeID') or ''
            tdt = txt(sub, 'TargetDataTypeID') or ''
            gid = txt(sub, 'LinkedBodyPartID') or ''
            gain = sub.find('Gain')
            ginfo = ''
            if gain is not None:
                ginfo = f"gain={gain.get('Type') or txt(gain,'ClassName')}"
                for p in ('A', 'B', 'C', 'D', 'MidPoint', 'Slope'):
                    v = actual(gain, p)
                    if v is not None:
                        ginfo += f" {p}={v}"
            subs.append(f"    sub {scls.split('.')[-1]}: DT={dt} TGT={tdt} body={gid} {ginfo}")
        print(f"  ADAPTER {name!r} ({cls.split('.')[-1]})")
        for s in subs:
            print(s)

# ---- stimuli: TonicStimulus nonzero neurons are already printed. The global Stimulus:
for el in root.iter('Stimulus'):
    print("\n=== GLOBAL STIMULUS ===", ET.tostring(el, encoding='unicode')[:700])
