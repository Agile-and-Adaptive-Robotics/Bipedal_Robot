"""Neural-class inventory of an AnimatLab .asim or .aproj (read-only).

Prints: neuron class histogram, synapse-type histogram, synapse type table,
per-neuron class sample (params), chart names + output files, SimEndTime.
Usage: D:\\Anaconda\\envs\\myo\\python.exe inventory.py <file> [--neurons]
"""
import io
import sys
import xml.etree.ElementTree as ET

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def txt(el, tag):
    x = el.find(tag)
    return x.text if x is not None else None


def act(el, tag):
    x = el.find(tag)
    if x is None:
        return None
    return x.get("Actual", x.get("Value"))


def main():
    path = sys.argv[1]
    show_neurons = "--neurons" in sys.argv
    root = ET.parse(path).getroot()
    is_aproj = root.tag == "Project"

    # --- neural modules (IntegrateFire) ---
    print(f"=== {path}  (root={root.tag}) ===")
    for mod in root.iter("Node"):
        cn = txt(mod, "ClassName") or ""
        if "NeuralModule" in cn and mod.find("SynapseTypes") is not None:
            print(f"MODULE class={cn}")
    for mod in root.iter():
        if mod.tag not in ("Node", "NeuralModule"):
            continue
        cn = txt(mod, "ClassName") or ""
        if "NeuralModule" not in cn or mod.find("SynapseTypes") is None:
            continue
        ts = act(mod, "TimeStep")
        print(f"  module timestep={ts}")

        # synapse types
        st_types = []
        for st in mod.find("SynapseTypes"):
            st_types.append(st)
        print(f"  synapse types: {len(st_types)}")
        hist = {}
        for st in st_types:
            cls = (txt(st, "ClassName") or "?").split(".")[-1]
            hist[cls] = hist.get(cls, 0) + 1
        for k, v in sorted(hist.items()):
            print(f"    {k}: {v}")
        for st in st_types:
            cls = (txt(st, "ClassName") or "?").split(".")[-1]
            name = txt(st, "Name") or txt(st, "Text") or "?"
            eq = act(st, "EquilibriumPotential")
            g0 = act(st, "SynapticConductance")
            amp = txt(st, "SynAmp")
            dec = act(st, "DecayRate")
            print(f"    TYPE id={txt(st,'ID')[:8]} cls={cls} name={name!r} SynAmp={amp} g0={g0} eq={eq} decay={dec}")

        # neurons
        neurons = []
        nodes_parent = mod.find("Nodes")
        if nodes_parent is None:
            nodes_parent = mod
        for nd in nodes_parent:
            cn = (txt(nd, "ClassName") or "")
            if "Neurons." not in cn:
                continue
            neurons.append(nd)
        hist = {}
        for nd in neurons:
            cls = (txt(nd, "ClassName") or "?").split(".")[-1]
            hist[cls] = hist.get(cls, 0) + 1
        print(f"  neurons: {len(neurons)}")
        for k, v in sorted(hist.items()):
            print(f"    {k}: {v}")
        if show_neurons:
            for nd in neurons:
                cls = (txt(nd, "ClassName") or "?").split(".")[-1]
                name = txt(nd, "Text") or txt(nd, "Name") or "?"
                print(
                    f"    N id={txt(nd,'ID')[:8]} cls={cls} name={name!r} "
                    f"rest={act(nd,'RestingPotential')} tc={act(nd,'TimeConstant')} "
                    f"thr={act(nd,'InitialThreshold')} gca={act(nd,'MaxCaConductance')} "
                    f"tonic={act(nd,'TonicStimulus')} ahp_g={act(nd,'AHP_Conductance')} "
                    f"ahp_tc={act(nd,'AHP_TimeConstant')} relacc={act(nd,'RelativeAccomodation')} "
                    f"acctc={act(nd,'AccomodationTimeConstant')}"
                )

    # --- connexions / links ---
    nconn = len(root.findall(".//Connexion"))
    nlink = len(root.findall(".//Link"))
    print(f"  connexions(.asim)={nconn} links(.aproj)={nlink}")

    # --- charts (asim) ---
    for ch in root.iter("DataChart"):
        name = txt(ch, "Name") or "?"
        print(f"  CHART {name!r} EndTime={act(ch,'EndTime')} collect={act(ch,'CollectTimeWindow')}")
    # --- SimEndTime ---
    for se in root.iter("SimEndTime"):
        print(f"  SimEndTime={se.get('Actual', se.get('Value', se.text))} raw={se.text}")
        break
    print()


if __name__ == "__main__":
    main()
