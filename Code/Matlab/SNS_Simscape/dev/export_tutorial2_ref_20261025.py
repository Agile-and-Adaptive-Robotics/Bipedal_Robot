# Export Tutorial-2 (5-neuron exc/inh/modulatory network) numpy reference.
import numpy as np
from scipy.io import savemat
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network

nrn = NonSpikingNeuron()
s_exc = NonSpikingSynapse(max_conductance=1.0, reversal_potential=50.0)
s_inh = NonSpikingSynapse(max_conductance=1.0, reversal_potential=-40.0)
s_mod = NonSpikingSynapse(max_conductance=1.0, reversal_potential=0.0)

net = Network(name='Tutorial 2 Network')
net.add_neuron(nrn, name='SourceNrn')
net.add_neuron(nrn, name='Dest1')
net.add_connection(s_exc, 'SourceNrn', 'Dest1')
net.add_neuron(nrn, name='Dest2')
net.add_neuron(nrn, name='Dest2In')
net.add_connection(s_exc, 'SourceNrn', 'Dest2')
net.add_connection(s_exc, 'Dest2', 'Dest2In')
net.add_connection(s_inh, 'Dest2In', 'Dest2')
net.add_neuron(nrn, name='Dest3')
net.add_connection(s_exc, 'SourceNrn', 'Dest3')
net.add_connection(s_mod, 'Dest1', 'Dest3')
net.add_input(dest='SourceNrn')
for n in ('SourceNrn', 'Dest1', 'Dest2', 'Dest2In', 'Dest3'):
    net.add_output(n)

dt, t_max = 0.01, 50.0
t = np.arange(0, t_max, dt)
inputs = np.zeros([len(t), 1]) + 20.0
model = net.compile(backend='numpy', dt=dt)
data = np.zeros([len(t), 5])
for i in range(len(t)):
    data[i, :] = model(inputs[i, :])
data = data.T
print("finals: Src %.4f Dest1 %.4f Dest2 %.4f Dest2In %.4f Dest3 %.4f mV"
      % tuple(data[:, -1]))
savemat(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\tutorial2_ref.mat",
        {"t": t, "Src": data[0], "Dest1": data[1], "Dest2": data[2],
         "Dest2In": data[3], "Dest3": data[4]}, do_compression=True)
print("saved results/tutorial2_ref.mat")
