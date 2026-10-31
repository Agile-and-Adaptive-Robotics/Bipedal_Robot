"""Independent conservative ground gate using recorded visual foot-mesh geometry."""
import numpy as np,mujoco
from scipy.spatial import ConvexHull

def check_recording(model,recording,samples,body_weight):
    data=mujoco.MjData(model);plane=model.geom('ground-plane').id;geometry=[]
    for side in ('r','l'):
        shapes=[]
        for g in range(model.ngeom):
            if model.body(int(model.geom_bodyid[g])).name not in ('calcn_'+side,'toes_'+side) or model.geom_type[g]!=mujoco.mjtGeom.mjGEOM_MESH:continue
            mesh=model.geom_dataid[g];a=model.mesh_vertadr[mesh];v=model.mesh_vert[a:a+model.mesh_vertnum[mesh]].copy()
            if len(v)>4:v=v[ConvexHull(v).vertices]
            shapes.append((g,v))
        if not shapes:raise ValueError('Foot clearance requires actual calcaneus/toes meshes')
        geometry.append(shapes)
    n=min(len(samples),len(recording['t']));idx=np.flatnonzero((recording['t'][:n]>=5)&(recording['t'][:n]<13))[::10];height=np.empty((len(idx),2))
    for row,i in enumerate(idx):
        data.qpos[:]=recording['qfull'][i];mujoco.mj_kinematics(model,data);normal=data.geom_xmat[plane].reshape(3,3)[:,2].copy()
        if normal[2]<0:normal=-normal
        for side,shapes in enumerate(geometry):height[row,side]=min(float(np.min(v@(data.geom_xmat[g].reshape(3,3).T@normal))+(data.geom_xpos[g]-data.geom_xpos[plane])@normal) for g,v in shapes)
    feet=[];times=recording['t'][idx]
    for side in (0,1):
        air=(height[:,side]>.001)&(samples[idx,7+side]/body_weight<.025);bounds=np.r_[0,np.flatnonzero(air[1:]!=air[:-1])+1,len(air)];episodes=[]
        for a,b in zip(bounds[:-1],bounds[1:]):
            if b-a>=4 and air[a] and height[a:b,side].max()>=.01:episodes.append(dict(start_s=float(times[a]),end_s=float(times[b-1]+.02),duration_s=float((b-a)*.02),peak_m=float(height[a:b,side].max())))
        feet.append(dict(sustained_swings=len(episodes),episodes=episodes,maximum_clearance_m=float(height[:,side].max()) if len(height) else None))
    return dict(accepted=all(f['sustained_swings']>=3 for f in feet),feet=feet,
                scope='Engineering clearance gate: 5-13s; >=3 swings/foot with force<.025BW, entire visual foot mesh >1mm for >=80ms, and peak>=10mm. Geometry sampled every20ms. Additional to unchanged movement/contact/frontal gates; not a biological standard.')
