from run_neural_robot import *
model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_abduction_diagnostic_v2.xml'))
for side in ('r','l'):
    data=mujoco.MjData(model);seed=np.zeros(3)
    for angle in np.linspace(-.85,.15,111):seed=knee_pose(model,data,side,float(angle),seed)
    atzero=knee_pose(model,data,side,0.,seed)
    data=mujoco.MjData(model);neutral=knee_pose(model,data,side,0.,np.zeros(3))
    print(side,'old continuation at neutral (deg)',np.rad2deg(atzero),'true neutral',np.rad2deg(neutral))
