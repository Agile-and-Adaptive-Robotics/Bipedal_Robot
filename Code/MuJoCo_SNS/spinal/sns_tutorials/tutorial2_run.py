# Auto-extracted from SNS-Toolbox tutorial notebook (cells in order).
# matplotlib .show() calls removed; render calls save to file instead.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def _save_current_fig(name):
    plt.gcf().savefig(r'D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\sns_tutorials' / name, dpi=130,
                      bbox_inches='tight')


# --- cell ---
# Add the library to the path
# If jupyter cannot find SNS-Toolbox
import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
# Import packages and modules for designing the network
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.networks import Network
from sns_toolbox.renderer import render

# Import packages and modules for simulating the network
import numpy as np
import matplotlib.pyplot as plt

# --- cell ---
# Define a non-spiking neuron and excitatory/inhibitory connections as in tutorial_1
neuron_type = NonSpikingNeuron()
synapse_excitatory = NonSpikingSynapse(max_conductance=1.0, reversal_potential=50.0)
synapse_inhibitory = NonSpikingSynapse(max_conductance=1.0, reversal_potential=-40.0)

# --- cell ---
# Define a modulatory synapse (relative reversal potential is 0)
synapse_modulatory = NonSpikingSynapse(max_conductance=1.0, reversal_potential=0.0)

# --- cell ---
net = Network(name='Tutorial 2 Network')

# --- cell ---
# Make a neuron to receive external input
net.add_neuron(neuron_type, name='SourceNrn', color='black')

# Make a neuron which is excited by the source neuron
net.add_neuron(neuron_type, name='Dest1', color='blue')
net.add_connection(synapse_excitatory, 'SourceNrn', 'Dest1')

render(net) # Display the network graph in this notebook

# --- cell ---
net.add_neuron(neuron_type, name='Dest2', color='orange')
net.add_neuron(neuron_type, name='Dest2In', color='green')
net.add_connection(synapse_excitatory, 'SourceNrn', 'Dest2')
net.add_connection(synapse_excitatory, 'Dest2', 'Dest2In')
net.add_connection(synapse_inhibitory, 'Dest2In', 'Dest2')

render(net)

# --- cell ---
net.add_neuron(neuron_type, name='Dest3', color='red')
net.add_connection(synapse_excitatory, 'SourceNrn', 'Dest3')
net.add_connection(synapse_modulatory, 'Dest1', 'Dest3')

render(net)

# --- cell ---
net.add_input(dest='SourceNrn', name='Input', color='white')

render(net)

# --- cell ---
net.add_output('SourceNrn', name='OutSourceNrn', color='grey')
net.add_output('Dest1', name='OutDest1', color='grey')
net.add_output('Dest2', name='OutDest2', color='grey')
net.add_output('Dest2In', name='OutDest2In', color='grey')
net.add_output('Dest3', name='OutDest3', color='grey')

render(net)

# --- cell ---
# Set simulation parameters
dt = 0.01 # Simulation timestep in ms
t_max = 50 # Max simulation time in ms

# Initialize a vector of timesteps
t = np.arange(0, t_max, dt)

# Initialize vectors which store the input to our network, and for data to be written to during simulation from outputs
inputs = np.zeros([len(t),1])+20.0  # Input vector must be 2d, even if second dimension is 1
data = np.zeros([len(t),5])

# --- cell ---
model = net.compile(backend='numpy', dt=dt)

# --- cell ---
for i in range(len(t)):
    data[i,:] = model(inputs[i,:])
data = data.transpose() # Rearrange the data axes for easier plotting

"""Plot the data"""
# First section
plt.figure()
plt.title('First Section')
plt.plot(t,data[:][0],label='SourceNrn',color='black')  # When plotting, all data needs to be transposed first
plt.plot(t,data[:][1],label='Dest1',color='blue')
plt.legend()

# Second section
plt.figure()
plt.title('Second Section')
plt.plot(t,data[:][0],label='SourceNrn',color='black')
plt.plot(t,data[:][2],label='Dest2',color='orange')
plt.plot(t,data[:][3],label='Dest2In',color='green')
plt.legend()

# Third section
plt.figure()
plt.title('Third Section')
plt.plot(t,data[:][0],label='SourceNrn',color='black')
plt.plot(t,data[:][1],label='Dest1',color='blue')
plt.plot(t,data[:][4],label='Dest3',color='red')
plt.legend()

# --- cell ---

