"""Probe: does sns_toolbox 1.5.2's numpy backend handle a MIXED network
(SpikingNeuron -> SpikingSynapse -> NonSpikingNeuron) without NaNs?
The compiler sets theta_0=float_max for non-spiking neurons but leaves
tau_theta=0 -> time_factor_threshold=inf -> inf*0=nan risk in the theta
update. Run: D:\\Anaconda\\envs\\myo\\python.exe probe_hybrid.py
"""
import numpy as np
from sns_toolbox.networks import Network
from sns_toolbox.neurons import SpikingNeuron, NonSpikingNeuron
from sns_toolbox.connections import SpikingSynapse

net = Network(name="hybrid probe")
net.add_neuron(SpikingNeuron(threshold_time_constant=5.0,
                             threshold_initial_value=1.0,
                             membrane_capacitance=0.005,
                             membrane_conductance=1.0,
                             resting_potential=0.0, bias=2.0),
               name="src")
net.add_neuron(NonSpikingNeuron(membrane_capacitance=0.03,
                                membrane_conductance=1.0,
                                resting_potential=0.0, bias=0.0),
               name="dst")
net.add_input("src", name="I0")
net.add_connection(SpikingSynapse(max_conductance=0.5,
                                  reversal_potential=8.0,
                                  time_constant=0.005,
                                  conductance_increment=0.05),
                   "src", "dst")
net.add_output("src", name="Ov", spiking=False)
net.add_output("src", name="Os", spiking=True)
net.add_output("dst", name="Od")

model = net.compile(backend="numpy", dt=0.0005)
print("spiking flag:", model.spiking)
print("time_factor_threshold:", model.time_factor_threshold)
print("theta_0:", model.theta_0)
n = 4000
out = np.zeros((n, 3))
nspk = 0
for k in range(n):
    o = model([0.0])
    out[k] = o
    nspk += int(model.spikes[0] == -1)
print("src spikes in 2 s:", nspk)
print("dst V finite:", np.all(np.isfinite(out[:, 2])),
      " min/max:", float(out[:, 2].min()), float(out[:, 2].max()))
print("model.V all finite:", np.all(np.isfinite(model.V)))
print("theta finite:", np.all(np.isfinite(model.theta)))
print("g_spike finite:", np.all(np.isfinite(model.g_spike)))
