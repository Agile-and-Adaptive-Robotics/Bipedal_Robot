"""Report reduced-coordinate torque capacity and tracking at a saved posture."""
import json,sys
from run_neural_robot import *
def main():
    model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'));data=mujoco.MjData(model)
    dynamic=DynamicRoutes(model,json.loads((HERE/'robotbody_routes.json').read_text()))
    trace=np.load(HERE/'runs/step_drive3.npz')
    for time_point in (0.,4.,18.):
        index=np.argmin(abs(trace['time']-time_point));data.qpos[:]=trace['qpos'][index];data.qvel[:]=trace['qvel'][index];data.act[:]=trace['activation'][index]
        dynamic.update(data);print('Time',time_point,'knees',dynamic.coordinates(data))
        for side in ('r','l'):
            j=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'thigh_'+side+'_joint');qa=model.jnt_qposadr[j];da=model.jnt_dofadr[j]
            print('Hip',side,'rotation',np.rad2deg(rotation_error(data.qpos[qa:qa+4])))
            for label,tangent in [('hipZ',np.eye(model.nv)[da+2])]:
                arms=dynamic.moment_arms[model.actuator_trnid[:,0]]@tangent;cap=-model.actuator_gainprm[:,0]*arms
                print(label,'gravity',data.qfrc_bias@tangent,'capacity',cap[cap<0].sum(),cap[cap>0].sum())
            joints=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,p+'_'+side+'_joint') for p in ('back','shank','front')]
            addresses=model.jnt_qposadr[joints];base=np.array(data.qpos);angle=dynamic.coordinates(data)['knee_angle_'+side]
            scratch=mujoco.MjData(model);scratch.qpos[:]=base
            plus=base.copy();plus[addresses]=knee_pose(model,scratch,side,angle+1e-5,base[addresses])
            minus=base.copy();minus[addresses]=knee_pose(model,scratch,side,angle-1e-5,base[addresses])
            tangent=np.zeros(model.nv);mujoco.mj_differentiatePos(model,tangent,2e-5,minus,plus)
            arms=dynamic.moment_arms[model.actuator_trnid[:,0]]@tangent;cap=-model.actuator_gainprm[:,0]*arms
            print('Knee',side,'gravity',data.qfrc_bias@tangent,'capacity',cap[cap<0].sum(),cap[cap>0].sum(),'actual muscle',data.qfrc_actuator@tangent)
        print('lumbar force cap',[(mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i),round(float(data.act[i]),3)) for i in range(model.nu) if 'ercspn' in mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i)])
if __name__=='__main__':main()
