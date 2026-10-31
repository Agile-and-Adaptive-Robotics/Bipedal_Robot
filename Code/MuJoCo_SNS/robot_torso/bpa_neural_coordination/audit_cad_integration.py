"""Inventory actual CAD route inputs; refuse unconfigured actuator simulation."""
import os,sys,json,hashlib
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
import mujoco
from cad_interface import REQUIRED,validate_configuration,RoutedBPABank
HERE=Path(__file__).resolve().parent
ROBOT=HERE.parent
sys.path.insert(0,str(ROBOT));from dynamic_routes import DynamicRoutes

def main():
    out=HERE/'runs'/'cad_parameter_audit_20261031';out.mkdir(parents=True,exist_ok=True)
    geometry_path=ROBOT/'robotbody_routes.json';model_path=ROBOT/'ah001_neural.xml'
    geometry=json.loads(geometry_path.read_text());model=mujoco.MjModel.from_xml_path(str(model_path));data=mujoco.MjData(model)
    dynamic=DynamicRoutes(model,geometry);dynamic.update(data)
    lookup={r['name']:r for r in geometry['routes']};rows=[];template={}
    names=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
    for i,name in enumerate(names):
        source=lookup[name];tid=int(model.actuator_trnid[i,0])
        assignment=source['actuator_assignment']
        rows.append(dict(name=name,current_route_length_m=float(data.ten_length[tid]),
            reference_capacity_N=float(-model.actuator_gainprm[i,0]),
            reference_capacity_provenance='Master MIF target, not measured BPA capacity',
            source_human_optimal_fiber_length_m=source['source_optimal_fiber_length_m'],
            source_human_tendon_slack_length_m=source['source_tendon_slack_length_m'],
            source_force_model=source['force_model'],master_name=assignment['master_name']))
        template[name]={key:None for key in REQUIRED}
        template[name]['reference_metadata_only']=rows[-1]
    manifest=dict(machine='easteregg2',status='UNCONFIGURED: fill route-specific dimensions and assumptions before use',
        command_mode='direct_pressure',actuators=template,
        notes=['Human optimal fiber lengths and tendon slack lengths are not BPA sleeve/cable dimensions.',
               'Master muscle MIFs are targets, not measured pneumatic capacity.',
               'Dynamic response parameters may be explicit assumptions, never silently called identified hardware values.',
               'Use the full moving-point length Jacobian and tensile virtual-work sign.'])
    missing=validate_configuration(manifest,names)
    assert len(missing)==54 and model.nu==54
    (HERE/'CAD_BPA_PARAMETER_TEMPLATE.json').write_text(json.dumps(manifest,indent=2))
    # Verify the proposed plumbing refuses this genuinely incomplete CAD config.
    refused=False
    try:RoutedBPABank(names,manifest,model.opt.timestep)
    except ValueError:refused=True
    assert refused
    # Tiny supported unit test of pressure state and virtual-work sign, no CAD stepping.
    fixture_entry=dict(diameter_mm=20,resting_length_m=.3,maximum_contraction_fraction=.25,nonmuscle_length_m=.06,
        series_stiffness_N_m=2e5,parallel_count=1,pressure_max_kpa=620.,tau_fill_s=.04,tau_vent_s=.06,delay_s=.004,
        provenance='Explicit synthetic fixture; not a CAD hardware assignment')
    bank=RoutedBPABank(['F','E'],dict(command_mode='direct_pressure',actuators={'F':fixture_entry,'E':fixture_entry}),.002)
    history=[]
    for k in range(100):history.append(bank.step([.5,.25],[.33,.33]).copy())
    assert bank.pressure[0]>bank.pressure[1]>0 and max(bank.pressure)<=620
    assert RoutedBPABank.generalized_force([[-.03],[.03]],[100.,0.])[0]==3.
    # Compare native MuJoCo tensile actuation with historical glue's positive sign.
    sign_model=mujoco.MjModel.from_xml_string('<mujoco><worldbody><body><joint name="j" axis="0 0 1"/><geom type="sphere" size=".1" mass="1"/></body></worldbody><tendon><fixed name="F"><joint joint="j" coef="-.03"/></fixed></tendon><actuator><general name="F" tendon="F" gainprm="-1"/></actuator></mujoco>')
    sd=mujoco.MjData(sign_model);sd.ctrl[0]=100.;mujoco.mj_forward(sign_model,sd)
    native=float(sd.qfrc_actuator[0]);historical=float((sd.actuator_moment.T@np.array([100.]))[0])
    assert native==3. and historical==-3.
    report=dict(model=str(model_path),source_routes=str(geometry_path),total_source_routes=len(geometry['routes']),
        active_actuators=model.nu,actuators_per_side={s:sum(n.endswith('_'+s) for n in names) for s in ('r','l')},
        mass_kg=float(model.body_mass.sum()),gravity=model.opt.gravity.tolist(),missing_parameter_routes=len(missing),
        required_fields=list(REQUIRED),route_inventory=rows,unconfigured_simulation_refused=refused,
        direct_pressure_bank_unit_check=True,virtual_work_sign_check=dict(native_tensile_generalized_force_Nm=native,
            positive_tension_times_length_gradient_Nm=historical,
            finding='Native virtual work requires -J_length^T*T. The pre-fix bpa_mujoco.py snapshot used the opposite sign; the helper is now corrected and covered by check_bpa_glue_parity.py. Historical results require explicit replay.'),
        source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (geometry_path,model_path,ROBOT.parent/'bpa_mujoco.py')})
    (out/'CAD_INTEGRATION_AUDIT.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k!='route_inventory'},indent=2))

if __name__=='__main__':main()
