"""Render Ben's BilateralRG connectome using the official SNS Toolbox diagram renderer.
v2: adds shape=[1] (API trap), fixes graphviz PATH."""
import json
import os
import sys
from pathlib import Path

# Graphviz binary lives in the myo env, not on system PATH
MYO_ENV = Path(r"C:\Users\Ben Bolen\.conda\envs\myo")
os.environ["PATH"] = str(MYO_ENV / "Library" / "bin") + os.pathsep + os.environ.get("PATH", "")

import numpy as np
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron, SpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse, SpikingSynapse

BEN_JSON = Path(r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models\Walker_2_Layer_CPG_BilateralRG_connectome.json")
OUT = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\slideshow\bilateralrg_sns_diagram")


def classify(label):
    l = str(label).lower()
    if "heel" in l or "toe" in l:
        return ("sensory", True)
    if "rg" in l and ("ext" in l or "flx" in l):
        return ("rg", False)
    if "pf" in l and ("hip" in l or "knee" in l or "ank" in l):
        return ("pf", False)
    if "mn" in l and "re" not in l.split()[-1].lower():
        return ("mn", False)
    if "mn" in l:
        return ("mn", False)  # Renshaw cells classified as MN-adjacent
    if " ia" in l or " ib" in l or " ii" in l or "_ia" in l or "_ib" in l or "_ii" in l:
        return ("sensory", False)
    if "ia in" in l or "ib in" in l:
        return ("in", False)
    if "in" in l.split() or "renshaw" in l or "_re" in l.lower():
        return ("in", False)
    if "stimulus" in l:
        return ("input", True)
    if "c1_" in l.lower() or "v3_" in l.lower():
        return ("commissural", False)
    return ("other", False)


def main():
    data = json.loads(BEN_JSON.read_text(encoding="utf8"))
    nodes = data.get("nodes", [])
    edges = data.get("edges", data.get("synapses", []))
    print(f"Ben's connectome: {len(nodes)} nodes, {len(edges)} edges")

    net = Network()
    ns_neuron = NonSpikingNeuron()
    spk_neuron = SpikingNeuron()

    import re
    label_to_pop = {}
    for n in nodes:
        label = str(n.get("label", n.get("name", n.get("id", "?"))))
        layer, is_spiking = classify(label)
        pop_name = re.sub(r"[^A-Za-z0-9_]", "_", label)[:30]
        neuron = spk_neuron if is_spiking else ns_neuron
        try:
            net.add_population(neuron, shape=[1], name=pop_name)
            label_to_pop[label] = pop_name
        except Exception as e:
            print(f"  SKIP: {label} ({type(e).__name__}: {e})")

    print(f"Created {len(label_to_pop)} populations")

    n_syn = 0
    for e in edges:
        src = str(e.get("from", e.get("source", "")))
        dst = str(e.get("to", e.get("target", "")))
        if src not in label_to_pop or dst not in label_to_pop:
            continue
        etype = str(e.get("type", e.get("sign", ""))).lower()
        is_inh = "inh" in etype
        g = float(e.get("conductance", e.get("weight", e.get("gain", 0.1))))
        try:
            syn = SpikingSynapse(g, is_inh)
            net.add_connection(syn, label_to_pop[src], label_to_pop[dst])
            n_syn += 1
        except Exception:
            try:
                syn = NonSpikingSynapse(g, is_inh)
                net.add_connection(syn, label_to_pop[src], label_to_pop[dst])
                n_syn += 1
            except Exception as ex2:
                print(f"  SKIP edge: {src}->{dst} ({ex2})")

    print(f"Created {n_syn} synapses")

    from sns_toolbox.renderer import render
    out_base = str(OUT.with_suffix(""))
    render(net, view=False, save=True, filename=out_base, img_format="png")
    print(f"Saved {out_base}.png")


if __name__ == "__main__":
    main()
