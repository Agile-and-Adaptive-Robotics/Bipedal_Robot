"""Runtime moving/conditional paths with complete virtual-work derivatives.

MuJoCo differentiates body motion with site_pos fixed. Moving-site terms are
added explicitly so muscle force and afferent velocity include guide motion.
Each tendon keeps its allocated path storage; only active site IDs and the
active count change. Other tendons' storage never moves.
"""
import numpy as np
import mujoco
from numba import njit


@njit(cache=True)
def spline_value_derivative(x,knots,coefficients,endpoints):
    if x<knots[0]:
        return endpoints[0,0]+endpoints[0,1]*(x-knots[0]),endpoints[0,1]
    if x>knots[-1]:
        return endpoints[1,0]+endpoints[1,1]*(x-knots[-1]),endpoints[1,1]
    index=min(np.searchsorted(knots,x,side='right')-1,len(coefficients)-1)
    index=max(0,index)
    t=x-knots[index];a,b,c,d=coefficients[index]
    return a+t*(b+t*(c+t*d)),b+t*(2*c+3*t*d)


class DynamicRoutes:
    def __init__(self,model,configuration):
        self.model=model
        self.routes=configuration['routes']
        self.coordinate_gradients={}
        self.qpos_terms={}
        for side in ('r','l'):
            terms=[];gradient=np.zeros(model.nv)
            for prefix in ('back','shank'):
                joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,prefix+'_'+side+'_joint')
                sign=float(model.jnt_axis[joint,2])
                if abs(abs(sign)-1)>1e-8: raise ValueError('Knee axis must follow OpenSim Z.')
                terms.append((int(model.jnt_qposadr[joint]),sign))
                gradient[model.jnt_dofadr[joint]]=sign
            self.qpos_terms['knee_angle_'+side]=terms
            self.coordinate_gradients['knee_angle_'+side]=gradient
        self.tendon_ids=[];self.site_ids=[];self.conditions=[];self.moving=[]
        self.active=[];self.addresses=[];self.capacities=[]
        for route in self.routes:
            tid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_TENDON,route['name'])
            if tid<0: raise ValueError('Missing tendon '+route['name'])
            self.tendon_ids.append(tid)
            ids=np.array([mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_SITE,p['id']) for p in route['points']],dtype=int)
            if np.any(ids<0): raise ValueError('Missing path site')
            self.site_ids.append(ids)
            self.addresses.append(int(model.tendon_adr[tid]))
            self.capacities.append(int(model.tendon_num[tid]))
            if len(ids)!=self.capacities[-1]: raise ValueError('Allocate all source path points before runtime switching.')
            conditional=[];moving=[]
            for index,point in enumerate(route['points']):
                if point.get('conditional_coordinate'):
                    name=point['conditional_coordinate']
                    if name not in self.qpos_terms: raise ValueError('Unsupported path coordinate '+name)
                    conditional.append((index,name,*point['conditional_range_rad']))
                functions=point.get('motion_functions',[])
                if functions:
                    records=[]
                    for f in functions:
                        if f['coordinate'] not in self.qpos_terms: raise ValueError('Unsupported moving coordinate '+f['coordinate'])
                        records.append((f['coordinate'],np.array(f['knots']),np.array(f['coefficients']),np.array(f['endpoints'])))
                    moving.append((index,ids[index],np.array(point['position_m']),np.array(point['source_location_parent_m']),
                                   np.array(point['motion_to_target_matrix']),records))
            self.conditions.append(conditional);self.moving.append(moving)
            self.active.append(np.ones(len(ids),dtype=bool))
        self.extra_length_jacobian=np.zeros((model.ntendon,model.nv))
        self.active_counts=np.array([len(x) for x in self.site_ids])
        self.motion_records=[]

    def coordinates(self,data):
        return {name:sum(data.qpos[address]*sign for address,sign in terms) for name,terms in self.qpos_terms.items()}

    def update(self,data):
        model=self.model;coordinates=self.coordinates(data)
        self.motion_records=[]
        for route_index,tid in enumerate(self.tendon_ids):
            mask=np.ones(len(self.site_ids[route_index]),dtype=bool)
            for point_index,name,lo,hi in self.conditions[route_index]:
                mask[point_index]=lo<=coordinates[name]<=hi
            ids=self.site_ids[route_index][mask]
            if len(ids)<2: raise ValueError('Conditional path has fewer than two points')
            address=self.addresses[route_index]
            model.wrap_objid[address:address+len(ids)]=ids
            model.tendon_num[tid]=len(ids)
            self.active[route_index]=mask;self.active_counts[route_index]=len(ids)
            for index,sid,base,source_base,mapping,functions in self.moving[route_index]:
                position=np.zeros(3);derivatives=[]
                for axis,(name,knots,coefficients,endpoints) in enumerate(functions):
                    position[axis],derivative=spline_value_derivative(coordinates[name],knots,coefficients,endpoints)
                    derivatives.append((axis,name,derivative))
                model.site_pos[sid]=base+mapping@(position-source_base)
                if mask[index]:
                    self.motion_records.append((route_index,tid,sid,mapping,derivatives))
        mujoco.mj_forward(model,data)
        self.extra_length_jacobian[:]=0
        for route_index,tid,sid,mapping,derivatives in self.motion_records:
            ids=self.site_ids[route_index][self.active[route_index]]
            points=data.site_xpos[ids]
            segments=points[1:]-points[:-1]
            norms=np.linalg.norm(segments,axis=1)
            units=segments/np.maximum(norms[:,None],1e-12)
            index=int(np.flatnonzero(ids==sid)[0]);gradient=np.zeros(3)
            if index>0: gradient+=units[index-1]
            if index<len(ids)-1: gradient-=units[index]
            body=model.site_bodyid[sid];rotation=data.xmat[body].reshape(3,3)
            for axis,name,derivative in derivatives:
                scalar=float(gradient@(rotation@mapping[:,axis]))*derivative
                self.extra_length_jacobian[tid]+=scalar*self.coordinate_gradients[name]
        self.length_jacobian=data.ten_J+self.extra_length_jacobian
        self.moment_arms=-self.length_jacobian
        self.velocities=self.length_jacobian@data.qvel
        return self

    def apply_force_correction(self,data):
        force=np.zeros(self.model.ntendon)
        np.add.at(force,self.model.actuator_trnid[:,0],data.actuator_force)
        data.qfrc_applied[:]+=self.extra_length_jacobian.T@force
