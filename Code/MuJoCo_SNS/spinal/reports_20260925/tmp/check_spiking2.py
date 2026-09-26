import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import numpy as np
from sns_toolbox.networks import Network
from sns_toolbox.neurons import SpikingNeuron
from sns_toolbox.connections import SpikingSynapse

net = Network(name="t")
net.add_neuron(SpikingNeuron(), name="s1")
net.add_neuron(SpikingNeuron(), name="s2")
net.add_connection(SpikingSynapse(max_conductance=2.0, reversal_potential=60.0), "s1", "s2")
net.add_input("s1")
net.add_output("s1", name="o1")
net.add_output("s2", name="o2")
compiled = net.compile(dt=0.0002, backend="numpy")
print("compiled type:", type(compiled).__name__)
o1, o2 = [], []
for i in range(25000):  # 5 s at 0.2 ms
    compiled.forward([10.0])
    o1.append(compiled.output[0])
    o2.append(compiled.output[1])
o1 = np.array(o1); o2 = np.array(o2)
print("s1 spikes:", int(np.sum(o1 > 0.5)), " s2 spikes (evoked):", int(np.sum(o2 > 0.5)))
print("first s1 spike t:", o1.argmax() * 0.0002)
print("first s2 spike t:", (o2 > 0.5).argmax() * 0.0002)
