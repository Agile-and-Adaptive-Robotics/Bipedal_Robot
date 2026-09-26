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

nm = root.find('.//NervousSystem/Node')
assert nm is not None and nm.find('Nodes') is not None

# module-level params + synapse types live on the IntegrateFire NeuralModule node
mod = None
for el in root.iter('Node'):
    if txt(el, 'ClassName') == 'IntegrateFireGUI.DataObjects.Behavior.NeuralModule' \
            and el.find('SynapseTypes') is not None:
        mod = el
        break
assert mod is not None
print("neural timestep (s):", actual(mod, 'TimeStep'),
      " AHPeq:", actual(mod, 'AHPEquilibriumPotential'),
      " CaEq:", actual(mod, 'CaEquilibriumPotential'),
      " SpikeStrength:", txt(mod, 'SpikeStrength'),
      " Refractory(ms):", actual(mod, 'RefractoryPeriod'))
print("neural timestep (s):", actual(nm, 'TimeStep'))

stypes = {}
for st in mod.find('SynapseTypes'):
    stypes[txt(st, 'ID')] = dict(
        cls=(txt(st, 'ClassName') or '').split('.')[-1],
        name=txt(st, 'Name'),
        eq=actual(st, 'EquilibriumPotential'),
        g0=actual(st, 'SynapticConductance'),
        decay=actual(st, 'DecayRate'),
    )
print("\n=== SYNAPSE TYPES ===")
for k, v in stypes.items():
    print(f"  {k[:8]}  {v['name']!r:28s} {v['cls']:18s} eq={v['eq']} g0={v['g0']} decay={v['decay']}")

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
    for gate in ('CaActivation', 'CaDeactivation'):
        g = nd.find(gate)
        if g is not None:
            pars[gate] = (actual(g, 'MidPoint'), actual(g, 'Slope'))
    neurons[nid] = dict(name=name, **pars)
    print(f"  {name!r:30s} {pars}")

print("\n=== CONNEXIONS (resolved) ===")
nbyid = {k: v['name'] for k, v in neurons.items()}
for lk in nm.find('Links'):
    o = txt(lk, 'OriginID')
    d = txt(lk, 'DestinationID')
    if o not in nbyid and d not in nbyid:
        continue
    st = txt(lk, 'SynapticTypeID')
    g = actual(lk, 'SynapticConductance')
    stv = stypes.get(st, {})
    print(f"  {nbyid.get(o,o)!r} -> {nbyid.get(d,d)!r}  [{stv.get('name')}] eq={stv.get('eq')} g={g}")

print("\n=== ADAPTERS ===")
for nd in nm.find('Nodes'):
    cls = txt(nd, 'ClassName') or ''
    if 'Adapter' not in cls:
        continue
    name = txt(nd, 'Text') or txt(nd, 'Name') or ''
    print(f"  ADAPTER {name!r} ({cls.split('.')[-1]})")
    for sub in (nd.find('Nodes') or []):
        scls = (txt(sub, 'ClassName') or '').split('.')[-1]
        dt = txt(sub, 'DataTypeID')
        tdt = txt(sub, 'TargetDataTypeID')
        gid = txt(sub, 'LinkedBodyPartID')
        gain = sub.find('Gain')
        ginfo = ''
        if gain is not None:
            gname = gain.get('Type') or txt(gain, 'ClassName') or txt(gain, 'Name')
            ginfo = f"gain={gname}"
            for p in ('A', 'B', 'C', 'D', 'MidPoint', 'Slope', 'Value'):
                v = actual(gain, p)
                if v is not None:
                    ginfo += f" {p}={v}"
        print(f"    sub {scls}: DT={dt} TGT={tdt} body={gid} {ginfo}")
