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
    
from sns_toolbox.networks import Network #, AdditionNetwork (This would import the code that we remake here
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingTransmissionSynapse
from sns_toolbox.networks import DivisionNetwork, MultiplicationNetwork, DifferentiatorNetwork
from sns_toolbox.networks import IntegratorNetwork
from sns_toolbox.renderer import render

# Let's define a custom functional subnetwork 'preset', in this case a network which takes a weighted sum of inputs

class AdditionNetwork(Network): # Note that this network is also a preset available from sns_toolbox.networks
    """
    Network which performs addition or subtraction of multiple inputs. Currently only supports non-spiking neurons.

    :param gains:       List of addition or subtraction weights.
    :type gains:        list, np.ndarray, or torch.tensor
    :param add_del_e:   Reversal potential of addition synapses, default is 100. Unit is millivolts (mV).
    :type add_del_e:    Number, optional
    :param sub_del_e:   Reversal potential of subtraction synapses, default is -40. Unit is millivolts (mV).
    :type sub_del_e:    Number, optional
    :param neuron_type: Neuron preset to use, default is sns_toolbox.design.neurons.NonSpikingNeuron.
    :type neuron_type:  sns_toolbox.design.neurons.NonSpikingNeuron, optional
    :param name:        Name of this network, default is 'Add'.
    :type name:         str, optional
    """
    def __init__(self,gains,add_del_e=100,sub_del_e=-40,neuron_type=NonSpikingNeuron(),name='Add', R=20.0, **kwargs):
        super().__init__(name=name,**kwargs)
        num_inputs = len(gains)
        self.add_neuron(neuron_type=neuron_type, name=name + 'Sum')
        for i in range(num_inputs):
            self.add_neuron(neuron_type, name=name + 'Src' + str(i))
            gain = gains[i]
            if gain > 0:
                conn = NonSpikingTransmissionSynapse(gain=gain, reversal_potential=add_del_e, e_lo=neuron_type.params['resting_potential'], e_hi=neuron_type.params['resting_potential']+R)
            else:
                conn = NonSpikingTransmissionSynapse(gain=gain, reversal_potential=sub_del_e, e_lo=neuron_type.params['resting_potential'], e_hi=neuron_type.params['resting_potential']+R)
            self.add_connection(conn, i + 1, name + 'Sum')

# Now let's import our network into another one, as we would normally use this functionality
net = Network(name='Tutorial 4 Network')

sum_net = AdditionNetwork([1,-1,-0.5,2])
net.add_network(sum_net, color='blue')
render(net)

# --- cell ---
# We can add more subnetworks
div_net = DivisionNetwork(1,0.5)
net.add_network(div_net, color='orange')
render(net)

# --- cell ---
mult_net = MultiplicationNetwork()
net.add_network(mult_net, color='green')
render(net)

# --- cell ---
diff_net = DifferentiatorNetwork()
net.add_network(diff_net,color='red')
render(net)

# --- cell ---
int_net = IntegratorNetwork()
net.add_network(int_net,color='purple')
render(net)

# --- cell ---

