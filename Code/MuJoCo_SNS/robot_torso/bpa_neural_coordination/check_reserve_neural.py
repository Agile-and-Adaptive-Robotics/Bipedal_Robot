"""Validate actual relay traces against SNS-Toolbox and actual MN-current formula."""
import json
from pathlib import Path
import numpy as np
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron
from reserve_neural_runner import DEST,F

def main():
 prefix=DEST/'w2lvar_fit_c0p05_h0p1'
 with np.load(prefix.with_suffix('.relay.npz')) as z:samples=z['samples'].copy();cap=float(z['current_cap_nA'])
 with np.load(prefix.with_suffix('.final_neural_state.npz')) as z:final=z['reserve_relay_voltage'].copy()
 net=Network(name='Explicit pressure reserve relays')
 for side in ('r','l'):
  name='RESERVE_pect_'+side;net.add_neuron(NonSpikingNeuron(membrane_capacitance=6.,membrane_conductance=1.,resting_potential=0.,bias=0.),name=name);net.add_input(name);net.add_output(name)
 model=net.compile(dt=2.,backend='numpy');maximum=0.
 for i,row in enumerate(samples):
  observed=model.forward(10.*row[:,3]);expected=samples[i+1,:,4] if i+1<len(samples) else final;maximum=max(maximum,float(np.max(abs(observed-expected))))
 assert maximum<1e-12,maximum
 with np.load(prefix.with_suffix('.motor_state.npz')) as z:u=z['excitation'].copy();names=list(z['actuator_names'])
 with np.load(prefix.with_suffix('.release_state.npz')) as z:initial=z['ctrl'].copy()
 ids=[names.index('pect_r'),names.index('pect_l')];old_mn=F.R.E_HI*np.vstack((initial,u[:-1]))[:,ids];rows=samples[300:300+len(u)]
 expected=np.clip(cap/40.*np.clip(rows[:,:,4]/10.,0,1)*(40.-old_mn),0.,cap);error=float(np.max(abs(expected-rows[:,:,5])));assert error<1e-12,error
 report=dict(machine='easteregg2',actual_source_trace=str(prefix),sns_toolbox_relay_max_error_mV=maximum,actual_homonymous_conductance_current_max_error_nA=error,nonzero_current_verified=bool(np.any(rows[:,:,5]>0)),current_cap_respected=bool(np.max(rows[:,:,5])<=cap),source_cell_count=2,biological_replication_claim=False)
 (DEST/'RELAY_UNIT_CHECKS.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
