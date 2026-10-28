from run_neural_robot import *
for filename in ('ah001_neural.xml','ah001_abduction_diagnostic.xml'):
    model=mujoco.MjModel.from_xml_path(str(HERE/filename));data=mujoco.MjData(model);mujoco.mj_forward(model,data)
    print(filename,'contacts',data.ncon)
    for c in data.contact:
        print(c.dist,[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,int(i)) for i in (c.geom1,c.geom2)])
for name in ('abd_feedback','abd_wide_heading8','abd_stop65'):
    trace=np.load(HERE/'runs'/ (name+'.npz'));print(name)
    for start in (0,2,5,15):
        mask=trace['time']>=start
        if mask.any():print(start,'knee min/max',np.rad2deg(trace['knee_coordinate'][mask].min(axis=0)),np.rad2deg(trace['knee_coordinate'][mask].max(axis=0)))
