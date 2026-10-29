"""SNS-Toolbox drawing of the separate four-cell PF candidate."""
import os,sys
from pathlib import Path
os.environ['PATH']=r'C:\Program Files\Graphviz\bin'+os.pathsep+os.environ.get('PATH','')
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron
from sns_toolbox.connections import NonSpikingSynapse
from sns_toolbox.renderer import render
H=Path(__file__).resolve().parent;O=H/'network_figures'
n=Network(name='Optional four-cell hip-knee PF module (per leg)')
for name in ('RF: hip flexion + knee extension','HAM: hip extension + knee flexion','IN from RF','IN from HAM'):n.add_neuron(NonSpikingNeuron(),name=name)
for source,target,g,rev in [('RF: hip flexion + knee extension','IN from RF',1.,8.),('IN from RF','HAM: hip extension + knee flexion',.3,-5.),('HAM: hip extension + knee flexion','IN from HAM',1.,8.),('IN from HAM','RF: hip flexion + knee extension',.3,-5.)]:
    n.add_connection(NonSpikingSynapse(max_conductance=g,reversal_potential=rev),source,target)
for fmt in ('svg','png'):render(n,view=False,save=True,filename=str(O/'biarticular_four_cell_candidate'),img_format=fmt)
