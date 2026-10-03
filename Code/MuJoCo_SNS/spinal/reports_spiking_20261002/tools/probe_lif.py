import numpy as np
from sns_toolbox.networks import Network
from sns_toolbox.neurons import SpikingNeuron
DT=0.0005
def mk(i, tt):
    n=Network()
    n.add_neuron(SpikingNeuron(threshold_time_constant=tt, threshold_initial_value=-50.0,
        threshold_proportionality_constant=0.0, threshold_leak_rate=1.0, threshold_increment=0.0,
        threshold_floor=-50.0, reset_potential=-60.0, membrane_capacitance=0.005,
        membrane_conductance=0.1, resting_potential=-70.0, bias=0.0), name="s")
    n.add_input("s")
    return n.compile(backend="numpy", dt=DT)
for tt in (100.0, 0.05):
    for I in (2.2, 2.5, 3.0, 5.0, 10.0, 21.9):
        m=mk(I,tt); sc=0
        for k in range(int(3.0/DT)):
            m([I])
            sc += int(m.spikes[0]==-1)
        print(f"tau_theta={tt} I={I:6.2f} f={sc/3.0:7.1f} Hz  V_last={m.V[0]:8.2f} theta={m.theta[0]:8.2f}")
