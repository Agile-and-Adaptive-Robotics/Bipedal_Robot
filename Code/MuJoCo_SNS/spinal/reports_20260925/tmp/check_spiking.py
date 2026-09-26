import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import sns_toolbox
print("sns_toolbox version:", getattr(sns_toolbox, "__version__", "?"))
from sns_toolbox import neurons
print("\nneurons module members:")
for nm in dir(neurons):
    if not nm.startswith("_"):
        print(" ", nm)
from sns_toolbox import connections
print("\nconnections module members:")
for nm in dir(connections):
    if not nm.startswith("_"):
        print(" ", nm)
# check numpy backend spiking support
from sns_toolbox.networks import Network
from sns_toolbox.neurons import SpikingNeuron
try:
    net = Network(name="t")
    net.add_neuron(SpikingNeuron(), name="s1")
    net.add_neuron(SpikingNeuron(), name="s2")
    from sns_toolbox.connections import SpikingSynapse
    net.add_connection(SpikingSynapse(max_conductance=2.0, reversal_potential=60.0), "s1", "s2")
    net.add_output("s2", name="out")
    compiled = net.compile(dt=0.0002, backend="numpy")
    import numpy as np
    spikes = 0
    v_out = []
    for i in range(5000):
        u = [10.0] + [0.0] * (len(net.get_inputs()) - 1) if hasattr(net, "get_inputs") else None
        compiled.forward([10.0, 0.0])
        v_out.append(compiled.output[0] if hasattr(compiled, "output") else 0)
    print("\nspiking numpy backend: OK, last outputs sample:", [round(float(x),3) for x in v_out[-5:]])
    print("n_spikes on s2 output:", sum(1 for v in v_out if float(v) > 0.5))
except Exception as e:
    print("\nspiking numpy backend FAILED:", type(e).__name__, e)
