"""Reproduce the SNS_Numpy failure at benchmark size 500 (spiking)."""
import numpy as np

from sns_toolbox.networks import Network
from sns_toolbox.neurons import SpikingNeuron
from sns_toolbox.connections import SpikingSynapse
from sns_toolbox import backends

net = Network()
s = SpikingNeuron(name="m", threshold_proportionality_constant=-1)
b = SpikingNeuron(name="bias", threshold_proportionality_constant=-1,
                  bias=10.0)
net.add_population(s, [42], name="ins")
net.add_population(b, [8], name="outs")
net.add_population(b, [6], name="connected")
net.add_population(b, [443], name="rest")
net.add_connection(SpikingSynapse(name="e"), "connected", "connected")
net.add_input("ins")
net.add_output("outs")
print("inputs:", net.get_num_inputs(), "populations:",
      net.get_num_populations())
try:
    m = backends.SNS_Numpy(net)
    print("compiled OK")
    m.forward(np.array([1.0]))
    print("stepped OK")
except Exception:
    import traceback
    traceback.print_exc()
