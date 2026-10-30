import os
from pathlib import Path
os.environ['PATH']=r'C:\Program Files\Graphviz\bin'+os.pathsep+os.environ.get('PATH','')
import graphviz
O=Path(__file__).resolve().parent/'network_figures'
g=graphviz.Digraph('balance',format='svg');g.attr(rankdir='LR',fontname='Arial',label='Actual balance dataflow: body-level command computation and SNS recruitment',labelloc='t')
g.attr('node',fontname='Arial',shape='box',style='rounded,filled',fillcolor='#FFF1D8')
g.node('com','Model CoM position / velocity\n(full-state simulation measurement)')
g.node('orientation','Torso orientation proxy\noptional weight; zero in selected PT case')
g.node('python','Python computation\nerror + filtering + 100 ms delay\nposition / velocity / optional acceleration gains')
g.node('force','Measured harness forces / torques\nrecorded for offline tuning')
g.node('tuning','Trial objective\nposture + loading + harness reaction')
g.node('stretch','Muscle length / stretch velocity\nlocal myostatic sensory signals',fillcolor='#E0EDF9')
g.node('spinal','Ia / II sensory pathways\nseparate from the PT branch',fillcolor='#E0EDF9')
g.node('sns','Four directional SNS task cells\nplus four inhibitory relays',shape='ellipse',fillcolor='#E0EDF9')
g.node('mn','Native motor-neuron populations',shape='ellipse',fillcolor='#E0EDF9')
g.edge('com','python');g.edge('orientation','python',style='dashed');g.edge('python','sns',label='computed input currents')
g.edge('sns','mn',label='fixed excitatory / inhibitory recruitment')
g.edge('stretch','spinal');g.edge('spinal','mn',label='enabled local reflex paths')
g.edge('force','tuning',style='dotted');g.edge('tuning','python',label='gain selection across trials',style='dotted')
g.render(str(O/'balance_dataflow'),cleanup=True)
