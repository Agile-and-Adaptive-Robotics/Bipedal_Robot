"""Export the tuned spinal network to JSON for the Simulink generator.

Composites the live config exactly like `runner --fitted --best`
(effective_tables("best")), mutates the params module IN PLACE (build_network
holds references to the dicts), builds the network with the same actuator
list the runner uses (all 92 model actuators — PRUNE_MUSCLES only zeroes
Fmax, it does not remove actuators), and walks the SNS Network object.

Output: spinal/spinal_net_export.json
  meta:   source note, dt, counts, constants (E_HI, reversal potentials)
  neurons:  name, tau_s, Vrest_mV, Gm_uS, Cm_nF (= 1000*tau)
  synapses: src, dst, g_uS, Esyn_mV, ThrPre_mV (0), SlopePre_mV (5)
  inputs:   port name -> destination neuron
  outputs:  MN neuron name -> actuator name (S(V) port feeds activation)
"""
from __future__ import annotations

import json
from pathlib import Path

import mujoco
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron

import build_network as bn
import params
from draw_circuit import effective_tables

HERE = Path(__file__).parent
MODEL = (HERE.parents[2] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
         / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3.xml")

# ---- composite the tuned configuration in place (runner --fitted --best) --
t = effective_tables("best")
params.G.clear(); params.G.update(t["G"])
params.TAU.clear(); params.TAU.update(t["TAU"])
for ph in params.W_PF_MN:
    params.W_PF_MN[ph].clear(); params.W_PF_MN[ph].update(t["W"][ph])
params.W_POSTURE.clear(); params.W_POSTURE.update(t["WPOST"])

# ---- build with the runner's actuator list -------------------------------
m = mujoco.MjModel.from_xml_path(str(MODEL))
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]
net = bn.build(acts, dt=params.DT, interleg=True)
n = net.net

# ---- dump -----------------------------------------------------------------
def neuronrec(pop):
    p = pop["type"].params
    assert isinstance(pop["type"], NonSpikingNeuron)
    tau = float(p["membrane_capacitance"]) / float(p["membrane_conductance"])
    return {
        "name": pop["name"],
        "tau_s": tau,
        "Vrest_mV": float(p["resting_potential"]),
        "Gm_uS": float(p["membrane_conductance"]),
        "Cm_nF": 1000.0 * float(p["membrane_capacitance"]),
        "bias_nA": float(p.get("bias", 0.0)),
    }

synapses = []
for c in n.connections:
    cp = dict(c["type"].params)
    cp.update(c["params"])            # instance overrides (e.g. scaled g)
    assert isinstance(c["type"], NonSpikingSynapse)
    synapses.append({
        "src": n.populations[c["source"]]["name"],
        "dst": n.populations[c["destination"]]["name"],
        "g_uS": float(cp["max_conductance"]),
        "Esyn_mV": float(cp["reversal_potential"]),
        "ThrPre_mV": float(cp.get("e_lo", 0.0)),
        "SlopePre_mV": float(cp.get("e_hi", 5.0)) - float(cp.get("e_lo", 0.0)),
    })

# The SNS Network stores every input port under the DEFAULT name 'Input';
# the real port names (the runner's u-vector semantics) are tracked in
# SpinalNetwork.inputs, appended in the same order as n.inputs.
assert len(n.inputs) == len(net.inputs), "input port count mismatch"
inputs = [{"port": net.inputs[i],
           "dst": n.populations[n.inputs[i]["destination"]]["name"]}
          for i in range(len(n.inputs))]

out = {
    "meta": {
        "source": t["note"],
        "dt_s": params.DT,
        "n_neurons": len(n.populations),
        "n_synapses": len(synapses),
        "n_inputs": len(inputs),
        "E_HI_mV": 5.0,
        "syn_thr_pre_mV": 0.0,
        "syn_slope_pre_mV": 5.0,
        "units": "mV, nA, uS, nF (SNS toolbox / SNS_Library conventions)",
    },
    "neurons": [neuronrec(p) for p in n.populations],
    "synapses": synapses,
    "inputs": inputs,
    "outputs": [{"mn": net.mn_names[a], "actuator": a, "actuator_id": i}
                for i, a in enumerate(acts)],
}
dest = HERE / "spinal_net_export.json"
dest.write_text(json.dumps(out, indent=1), encoding="utf-8")

# ---- sanity ---------------------------------------------------------------
idx = net.idx
assert out["meta"]["n_neurons"] == len(idx), "population/index mismatch"
for s in synapses:
    assert s["src"] in idx and s["dst"] in idx
exc = sum(1 for s in synapses if s["Esyn_mV"] > 0)
taus = sorted({round(r["tau_s"], 4) for r in out["neurons"]})
print(f"wrote {dest.name}: {out['meta']['n_neurons']} neurons, "
      f"{len(synapses)} synapses ({exc} exc / {len(synapses) - exc} inh), "
      f"{len(inputs)} inputs, {len(out['outputs'])} MN->actuator outputs")
print(f"source: {t['note']}")
print(f"distinct taus (s): {taus}")
