"""Independent source-SNS leak, conductance-current, and complete core-state checks."""
import json
from types import SimpleNamespace
import numpy as np
import mn_recruitment_experiment as E
from capacity_estimator_checks import Inner,context

def main():
 E.GAIN=8.;r=E.Recruitment(Inner(),['MN_pect_r','MN_pect_l']);c=context();c['adapter'].bank.scales=np.array([100.,100.]);raw=np.array([[60.,10.,.1751,.001],[65.,12.,.1751,.001]]);r.outer=SimpleNamespace(context=c,initialized=True,buffer=np.tile(raw,(4,1,1)),buffer_count=4,delay_steps=3,previous_force=np.array([10.,12.]),routes=['pect_r','pect_l'])
 from sns_toolbox.networks import Network
 from sns_toolbox.neurons import NonSpikingNeuron
 n=Network(name='Independent recruitment relay')
 for name in ('r','l'):n.add_neuron(NonSpikingNeuron(membrane_capacitance=6.,membrane_conductance=1.),name=name);n.add_input(dest=name);n.add_output(source=name)
 model=n.compile(dt=2.,backend='numpy');ve=0.;ie=0.;x=np.array([.1,.2])
 for step in range(96):
  vm=r.inner.c.V.copy();r.forward(x);row=r.logs[-1];expected=np.clip((8./40.)*np.clip(row[:,5]/10.,0,1)*(40.-vm),0.,8.);ie=max(ie,float(np.max(np.abs(expected-row[:,6]))));model.forward(10.*row[:,4]);ve=max(ve,float(np.max(np.abs(model.V-r.relay))))
 state=r.export_state()
 for _ in range(32):r.forward(x)
 final=r.export_state();r.import_state(state)
 for _ in range(32):r.forward(x)
 repeated=r.export_state();checks={k:bool(np.array_equal(final[k],repeated[k])) for k in final};assert ve==0. and ie==0. and all(checks.values())
 result=dict(independent_SNS_voltage_max_error_mV=ve,conductance_current_max_error_nA=ie,full_core_repeat=checks,signal_uses_no_actual_gain=True,unit_scope='Synthetic sensor/circuit fixture only; no standing inference. Whole transformed native model restoration additionally checked automatically during real recruitment.')
 (E.DEST/'INDEPENDENT_RECRUITMENT_CHECKS.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(SNS_error_mV=ve,current_error_nA=ie,all_states_exact=all(checks.values()))))

if __name__=='__main__':main()
