"""Independent SNS relay and explicit pressure-output-map validation on actual data."""
import json,argparse
from pathlib import Path
import numpy as np
from numba import njit
from actuator_interface import pressure_step
from sns_toolbox.networks import Network
from sns_toolbox.neurons import NonSpikingNeuron
from reserve_valve_runner import DEST

@njit(cache=True)
def hydraulic_replay(commands,pressure,history,count,delay,fill,vent):
 out=np.empty_like(commands)
 for k in range(len(commands)):
  for i in range(len(pressure)):
   history[i,count%history.shape[1]]=commands[k,i]
   reference=history[i,(count-delay[i])%history.shape[1]] if count>=delay[i] else 0.
   pressure[i]=pressure_step(pressure[i],reference,.002,fill[i],vent[i]);out[k,i]=pressure[i]
  count+=1
 return out

def main(tag='fit_p0p5_r0'):
 prefix=DEST/('w2lvar_'+tag)
 with np.load(prefix.with_suffix('.relay.npz')) as z:rows=z['samples'];cap=float(z['pressure_cap_kpa'])
 with np.load(prefix.with_suffix('.final_neural_state.npz')) as z:final=z['reserve_relay_voltage']
 net=Network(name='Explicit bilateral reserve-to-valve relays')
 for side in ('r','l'):
  name='RESERVE_pect_'+side;net.add_neuron(NonSpikingNeuron(membrane_capacitance=6.,membrane_conductance=1.,resting_potential=0.,bias=0.),name=name);net.add_input(name);net.add_output(name)
 model=net.compile(dt=2.,backend='numpy');error=0.
 for i,row in enumerate(rows):
  observed=model.forward(10.*row[:,3]);expected=rows[i+1,:,4] if i+1<len(rows) else final;error=max(error,float(np.max(abs(observed-expected))))
 assert error<1e-12 and np.max(10*rows[:,:,3])<=10. and np.min(rows[:,:,3])>=0.
 assert np.count_nonzero(rows[:,:,5])==0,'Valve path must not inject existing MN currents'
 with np.load(prefix.with_suffix('.valve.npz')) as z:valve=z['samples']
 expected=cap*np.clip(rows[300:300+len(valve),:,4]/10.,0,1);map_error=float(np.max(abs(expected-valve[:,:,1])))
 assert map_error==0. and np.max(valve[:,:,2])<=cap and np.min(valve[:,:,2])>=0.
 assert np.array_equal(valve[:,:,3],valve[:,:,0]+valve[:,:,2])
 r=json.loads(prefix.with_suffix('.result.json').read_text());cfg=json.loads(Path(r['bpa_configuration']).read_text());names=list(cfg['actuators']);fill=np.array([cfg['actuators'][a]['tau_fill_s'] for a in names]);vent=np.array([cfg['actuators'][a]['tau_vent_s'] for a in names]);delay=np.array([round(cfg['actuators'][a]['delay_s']/.002) for a in names],dtype=np.int64)
 with np.load(prefix.with_suffix('.release_state.npz')) as z:initial=z['pressure'].copy();history=z['pressure_history'].copy();count=int(z['pressure_step_count'])
 with np.load(prefix.with_suffix('.pressure.npz')) as z:command=z['samples'][:,0,:];actual=z['samples'][:,1,:]
 replay=hydraulic_replay(command,initial,history,count,delay,fill,vent);assert np.array_equal(replay,actual),'Actual pressure deviates from original delayed first-order physics'
 report=dict(machine='easteregg2',sns_toolbox_voltage_error_mV=error,voltage_to_pressure_output_error_kpa=map_error,existing_MN_current_injection_zero=True,relay_input_current_cap_nA=10.,output_increment_cap_kpa=cap,actual_pressure_state_not_modified=True,full86route_pressure_replay_exact=True,actual_pressure_replay_max_error_kpa=float(np.max(abs(replay-actual))),source=str(prefix))
 (DEST/'VALVE_NEURAL_UNIT_CHECKS.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--tag',default='fit_p0p5_r0');main(p.parse_args().tag)
