"""Independently review optional saved current accounting; never advance a model."""
from pathlib import Path
import argparse,json
import numpy as np


def review(source):
    source=Path(source);z=np.load(source/'aligned_trace.npz',allow_pickle=False)
    meta=json.loads((source/'ALIGNED_TRACE_METADATA.json').read_text(encoding='utf-8'))
    result={'source':str(source.resolve()),'samples':len(z['t_state']),
            'neural_cells':z['neural_names'].tolist(),'neural_and_state_times_exact':bool(np.array_equal(z['t_state'],z['t_neural'])),
            'all_neural_voltages_finite':bool(np.isfinite(z['neural_voltage']).all()),
            'maximum_regional_sum_error_N':float(abs(z['regional_vertical_N'].sum(axis=2)-z['foot_force_anatomical_N'][:,:,1]).max()),
            'current_metadata':meta.get('current_telemetry'),
            'scope':'Read-only audit of saved samples. No reconstruction of missing current, neuronal step, or physics operation.'}
    if 'neural_current_available'not in z:
        result['optional_currents_recorded']=False;return result
    available=z['neural_current_available'];codes,counts=np.unique(z['neural_current_status'],return_counts=True)
    result.update(optional_currents_recorded=True,available_samples=int(available.sum()),
                  all_samples_current_available=bool(available.all()),status_counts={str(int(k)):int(v)for k,v in zip(codes,counts)})
    if not available.any():return result
    voltage=z['neural_voltage'][available];last=z['neural_V_last_mV'][available]
    total=sum(z['neural_'+part+'_current_nA'][available]for part in('leak','bias','chemical','applied','ionic'))
    reconstructed=last+z['neural_time_factor_membrane'][available]*total
    residual=voltage-reconstructed
    result.update(independent_maximum_update_error_mV=float(abs(residual).max()),
                  saved_maximum_update_residual_mV=float(abs(z['neural_update_residual_mV'][available]).max()),
                  independent_closure_passed=bool(np.isfinite(residual).all()and abs(residual).max()<=1e-11))
    result['cells']={}
    for i,name in enumerate(z['neural_names'].tolist()):
        fields={}
        for part in('applied','chemical','leak','bias','ionic'):
            values=z['neural_'+part+'_current_nA'][available,i]
            fields[part+'_nA']={'minimum':float(values.min()),'median':float(np.median(values)),'maximum':float(values.max())}
        fields['maximum_update_error_mV']=float(abs(residual[:,i]).max());result['cells'][name]=fields
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',required=True,type=Path);a=p.parse_args()
    result=review(a.source);out=a.source/'SPARSE_CURRENT_REVIEW.json';out.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items()if k not in('cells','neural_cells','current_metadata')},indent=2))
if __name__=='__main__':main()
