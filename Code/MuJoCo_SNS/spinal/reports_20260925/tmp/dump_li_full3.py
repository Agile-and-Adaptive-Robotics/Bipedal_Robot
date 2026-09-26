import io, sys
import xml.etree.ElementTree as ET
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

path = r'D:\Github\Bipedal_Robot\Neuromechanical_Models\Li Model\walk tester rearranged.aproj'
root = ET.parse(path).getroot()

def txt(el, tag):
    x = el.find(tag)
    return x.text if x is not None else None

def act(el, tag):
    """Actual value of a triplet stored as ATTRIBUTES on the child element."""
    x = el.find(tag)
    if x is None:
        return None
    return x.get('Actual', x.get('Value'))

nm = root.find('.//NervousSystem/Node')
mod = None
for el in root.iter('Node'):
    if txt(el, 'ClassName') == 'IntegrateFireGUI.DataObjects.Behavior.NeuralModule' \
            and el.find('SynapseTypes') is not None:
        mod = el
        break
print("MODULE timestep(s):", act(mod, 'TimeStep'),
      " AHPeq:", act(mod, 'AHPEquilibriumPotential'),
      " CaEq:", act(mod, 'CaEquilibriumPotential'),
      " SpikeStrength:", txt(mod, 'SpikeStrength'),
      " SpikePeak:", act(mod, 'SpikePeak'),
      " Refractory(s):", act(mod, 'RefractoryPeriod'))

stypes = {}
for st in mod.find('SynapseTypes'):
    stypes[txt(st, 'ID')] = dict(
        cls=(txt(st, 'ClassName') or '').split('.')[-1], name=txt(st, 'Name'),
        eq=act(st, 'EquilibriumPotential'), g0=act(st, 'SynapticConductance'),
        decay=act(st, 'DecayRate'))
print("\n=== SYNAPSE TYPES ===")
for k, v in stypes.items():
    print(f"  {k[:8]}  {v['name']!r:32s} {v['cls']:18s} eq={v['eq']} g0={v['g0']} decay={v['decay']}")

print("\n=== NEURONS ===")
neurons = {}
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if 'Neurons.' not in cls:
        continue
    nid = txt(nd, 'ID')
    name = txt(nd, 'Text') or txt(nd, 'Name') or ''
    pars = dict(cls=cls.split('.')[-1], rest=act(nd, 'RestingPotential'),
                tc=act(nd, 'TimeConstant'), thr=act(nd, 'InitialThreshold'),
                acc=act(nd, 'RelativeAccomodation'), acctc=act(nd, 'AccomodationTimeConstant'),
                ahp_g=act(nd, 'AHP_Conductance'), ahp_tc=act(nd, 'AHP_TimeConstant'),
                ca_g=act(nd, 'MaxCaConductance'), tonic=act(nd, 'TonicStimulus'))
    for gate in ('CaActivation', 'CaDeactivation'):
        g = nd.find(gate)
        if g is not None:
            pars[gate] = (act(g, 'MidPoint'), act(g, 'Slope'))
    neurons[nid] = dict(name=name, **pars)
    print(f"  {name!r:30s} {pars}")

print("\n=== CONNEXIONS (resolved, with per-link g) ===")
nbyid = {k: v['name'] for k, v in neurons.items()}
n_conn = 0
for lk in nm.find('Links'):
    o, d = txt(lk, 'OriginID'), txt(lk, 'DestinationID')
    if o not in nbyid and d not in nbyid:
        continue
    st = txt(lk, 'SynapticTypeID')
    g = act(lk, 'SynapticConductance')
    stv = stypes.get(st, {})
    n_conn += 1
    print(f"  {nbyid.get(o,o)!r} -> {nbyid.get(d,d)!r}  [{stv.get('name')}] g={g}")

print("\n=== ADAPTERS (with subnode gains) ===")
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if 'Adapter' not in cls:
        continue
    name = txt(nd, 'Text') or txt(nd, 'Name') or ''
    print(f"  ADAPTER {name!r} ({cls.split('.')[-1]}) id={txt(nd,'ID')}")
    subs = nd.find('Nodes')
    if subs is None:
        print("    (no subnodes)")
        continue
    for sub in subs:
        scls = (txt(sub, 'ClassName') or '').split('.')[-1]
        dt, tdt, gid = txt(sub, 'DataTypeID'), txt(sub, 'TargetDataTypeID'), txt(sub, 'LinkedBodyPartID')
        gain = sub.find('Gain')
        ginfo = ''
        if gain is not None:
            gname = gain.get('Type') or txt(gain, 'Name') or txt(gain, 'ClassName')
            ginfo = f"gain={gname}"
            for p in ('A', 'B', 'C', 'D', 'MidPoint', 'Slope', 'Value'):
                v = act(gain, p)
                if v is not None:
                    ginfo += f" {p}={v}"
        print(f"    sub {scls}: DT={dt} TGT={tdt} body={gid} {ginfo}")

print("\n=== CONTACT SENSOR bodies (RigidBody nodes w/ ContactSensor) ===")
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if cls.endswith('.RigidBody'):
        print(f"  {txt(nd,'Text')!r} id={txt(nd,'ID')} partid={txt(nd,'LinkedBodyPartID')}")
        for sub in (nd.find('Nodes') or []):
            scls = (txt(sub, 'ClassName') or '').split('.')[-1]
            dt = txt(sub, 'DataTypeID')
            gain = sub.find('Gain')
            ginfo = ''
            if gain is not None:
                gname = gain.get('Type') or txt(gain, 'Name') or txt(gain, 'ClassName')
                ginfo = f"gain={gname}"
                for p in ('A', 'B', 'C', 'D', 'MidPoint', 'Slope'):
                    v = act(gain, p)
                    if v is not None:
                        ginfo += f" {p}={v}"
            print(f"    sub {scls} DT={dt} {ginfo}")

print("\n=== JOINT nodes ===")
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if cls.endswith('.Joint'):
        print(f"  {txt(nd,'Text')!r} id={txt(nd,'ID')} partid={txt(nd,'LinkedBodyPartID')}")
        for sub in (nd.find('Nodes') or []):
            scls = (txt(sub, 'ClassName') or '').split('.')[-1]
            dt = txt(sub, 'DataTypeID')
            gain = sub.find('Gain')
            ginfo = ''
            if gain is not None:
                gname = gain.get('Type') or txt(gain, 'Name') or txt(gain, 'ClassName')
                ginfo = f"gain={gname}"
                for p in ('A', 'B', 'C', 'D', 'MidPoint', 'Slope'):
                    v = act(gain, p)
                    if v is not None:
                        ginfo += f" {p}={v}"
            print(f"    sub {scls} DT={dt} {ginfo}")

print("\n=== MUSCLE nodes (behavior) ===")
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if cls.endswith('.Muscle'):
        print(f"  {txt(nd,'Text')!r} id={txt(nd,'ID')} partid={txt(nd,'LinkedBodyPartID')}")
