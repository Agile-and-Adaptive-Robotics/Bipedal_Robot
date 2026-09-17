"""Probe the compiled backend object's actual attributes."""
import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.networks import Network

n = Network()
n.add_neuron(NonSpikingNeuron(membrane_capacitance=0.05,
                              membrane_conductance=1.0,
                              resting_potential=0.0, bias=0.0), name="A")
n.add_input("A")
net = n.compile(backend="numpy", dt=0.002)
print(type(net))
print("module:", type(net).__module__)
attrs = [a for a in dir(net) if not a.startswith("__")]
print([a for a in attrs if a in ("u", "V", "num_inputs", "forward",
                                 "simulate", "step")])
print("has u:", hasattr(net, "u"), "| has V:", hasattr(net, "V"))
import inspect
print(inspect.getsourcefile(type(net)))
