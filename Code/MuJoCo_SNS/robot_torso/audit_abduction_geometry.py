"""Check recorded inter-leg contacts, knee stops, and equality closure."""
import argparse,json
from run_neural_robot import *
def audit(path):
    summary=json.loads(path.with_suffix('.json').read_text());trace=np.load(path)
    model=mujoco.MjModel.from_xml_path(summary.get('model_path',str(HERE/'ah001_neural.xml')));data=mujoco.MjData(model)
    minimum=0.;contacts=0;maximum_equality=0.
    for q in trace['qpos'][trace['time']>=15][::5]:
        data.qpos[:]=q;mujoco.mj_forward(model,data)
        mask=data.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY)
        if mask.any():maximum_equality=max(maximum_equality,float(np.abs(data.efc_pos[mask]).max()))
        for contact in data.contact:
            names=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_BODY,int(model.geom_bodyid[g])) or '' for g in (contact.geom1,contact.geom2)]
            if names[0].endswith(('_r','_l')) and names[1].endswith(('_r','_l')) and names[0][-1]!=names[1][-1]:
                contacts+=1;minimum=min(minimum,float(contact.dist))
    window=trace['time']>=15
    result=dict(run=path.stem,interleg_contact_samples=contacts,maximum_interleg_penetration_m=-minimum,
        maximum_knee_loop_position_error_m=maximum_equality,
        post_ramp_knee_min_deg=np.rad2deg(trace['knee_coordinate'][window].min(axis=0)).tolist(),
        post_ramp_knee_max_deg=np.rad2deg(trace['knee_coordinate'][window].max(axis=0)).tolist(),
        scope='Opposite-leg CAD mesh contacts and explicit relative-knee stop surrogates; no hardware-clearance validation')
    path.with_name(path.stem+'_geometry_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path)
    audit(p.parse_args().run)
