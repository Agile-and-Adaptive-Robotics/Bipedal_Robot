"""Hunt C/D hypothesis: joint primitives encode net active task torque, not coactivation."""
import os,runpy
from pathlib import Path
os.environ['HUNT_CAMPAIGN']='net_torque_pf_20261030_v1';os.environ['HUNT_PF_TARGET']='net_active_torque'
runpy.run_path(str(Path(__file__).with_name('train_fourcell_multicondition.py')),run_name='__main__')
