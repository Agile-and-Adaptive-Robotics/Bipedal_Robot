"""Extract the existing symmetric robot-body OpenSim routes without edits."""
import hashlib
import json
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
import copy
import numpy as np

os.environ['CONDA_PREFIX'] = str(Path(sys.executable).parent)
import opensim as osim

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
SOURCE = REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim'

def export_function(function,xml_node,coordinate):
    xnode=xml_node.find('.//SimmSpline/x')
    if xnode is None:
        raise ValueError('Expected source SimmSpline; do not approximate an unknown function.')
    knots=np.fromstring(xnode.text,sep=' ')
    vector=osim.Vector(1,0.)
    derivative=osim.StdVectorInt(); derivative.append(0)
    def evaluate(x):
        vector.set(0,float(x))
        return float(function.calcValue(vector)),float(function.calcDerivative(derivative,vector))
    coefficients=[]; error=0.
    for a,b in zip(knots[:-1],knots[1:]):
        ya,da=evaluate(a);yb,db=evaluate(b);h=b-a
        c=[ya,da,(3*(yb-ya)/h-2*da-db)/h,(2*(ya-yb)/h+da+db)/(h*h)]
        coefficients.append(c)
        for t in (.17,.53,.89):
            dx=t*h
            error=max(error,abs(evaluate(a+dx)[0]-sum(c[i]*dx**i for i in range(4))))
    assert error<1e-10,error
    endpoints=[evaluate(knots[0]),evaluate(knots[-1])]
    # Verify OpenSim's linear endpoint extrapolation independently.
    for index,direction in ((0,-1),(-1,1)):
        x=knots[index]+direction*.01
        y,d=endpoints[0 if index==0 else 1]
        assert abs(evaluate(x)[0]-(y+d*(x-knots[index])))<1e-9
    return dict(coordinate=coordinate.split('/')[-1],knots=knots.tolist(),coefficients=coefficients,
                endpoints=endpoints,maximum_value_check_error_m=error)


def vector(v):
    return [float(v.get(i)) for i in range(3)]


def main():
    model = osim.Model(str(SOURCE))
    model.setUseVisualizer(False)
    state = model.initSystem()
    coords = model.getCoordinateSet()
    for i in range(coords.getSize()):
        coords.get(i).setValue(state, 0., False)
    model.realizePosition(state)
    xml = ET.parse(SOURCE).getroot()
    xml_muscles = {x.get('name'):x for x in xml.iter('Thelen2003Muscle')}
    output = dict(source=str(SOURCE.relative_to(REPO)), sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                  source_status='Direct gait2392_robotbody right-side extraction; left paths explicitly mirrored from right in this inventory',
                  opensim_version=osim.GetVersionAndDate(), coordinates={coords.get(i).getName():coords.get(i).getValue(state) for i in range(coords.getSize())},
                  bodies={}, joints={}, muscles=[])
    bodies = model.getBodySet()
    for i in range(bodies.getSize()):
        body = bodies.get(i)
        transform = body.getTransformInGround(state)
        rot=transform.R()
        output['bodies'][body.getName()] = dict(origin_ground_m=vector(transform.p()),
            rotation_to_ground=[[float(rot.get(a,b)) for b in range(3)] for a in range(3)])
    joints = model.getJointSet()
    for i in range(joints.getSize()):
        joint = joints.get(i)
        output['joints'][joint.getName()] = dict(parent_frame_origin_ground_m=vector(joint.getParentFrame().getTransformInGround(state).p()),
            child_frame_origin_ground_m=vector(joint.getChildFrame().getTransformInGround(state).p()))
    muscles = model.getMuscles()
    for i in range(muscles.getSize()):
        muscle=muscles.get(i)
        name=muscle.getName()
        if not name.endswith('_r'):
            continue
        path=muscle.getGeometryPath()
        points=path.getPathPointSet()
        xml_points=list(xml_muscles[name].find('GeometryPath/PathPointSet/objects'))
        item=dict(name=name, max_isometric_force_N=float(muscle.getMaxIsometricForce()),
                  optimal_fiber_length_m=float(muscle.getOptimalFiberLength()),
                  tendon_slack_length_m=float(muscle.getTendonSlackLength()),
                  source_neutral_path_length_m=float(path.getLength(state)), points=[],
                  source_neutral_moment_arms_m={coord:float(muscle.computeMomentArm(state,coords.get(coord)))
                    for coord in ('hip_flexion_r','hip_adduction_r','hip_rotation_r','knee_angle_r',
                                  'hip_flexion_l','hip_adduction_l','hip_rotation_l','knee_angle_l','lumbar_extension')})
        for j in range(points.getSize()):
            point=points.get(j)
            parent=point.getParentFrame()
            source_xml=xml_points[j]
            assert source_xml.get('name')==point.getName()
            condition=source_xml.findtext('socket_coordinate')
            bounds=source_xml.findtext('range')
            enabled=True
            if source_xml.tag=='ConditionalPathPoint':
                coordinate=condition.split('/')[-1]
                lo,hi=map(float,bounds.split())
                enabled=lo <= coords.get(coordinate).getValue(state) <= hi
            local=point.getLocation(state)
            item['points'].append(dict(name=point.getName(), type=point.getConcreteClassName(),
                parent_frame=str(parent.getAbsolutePathString()), body=parent.findBaseFrame().getName(),
                location_parent_m=vector(local), location_ground_m=vector(parent.findStationLocationInGround(state,local)),
                neutral_enabled=enabled, condition_coordinate=condition, condition_range=bounds,
                source_xml=ET.tostring(source_xml,encoding='unicode')))
            if source_xml.tag=='MovingPathPoint':
                moving=osim.MovingPathPoint.safeDownCast(point)
                item['points'][-1]['motion_functions']=[export_function(getattr(moving,'get_'+axis+'_location')(),
                    source_xml.find(axis+'_location'),source_xml.findtext('socket_'+axis+'_coordinate')) for axis in 'xyz']
        wrapset=path.getWrapSet()
        item['path_wrap_count']=wrapset.getSize()
        output['muscles'].append(item)
    for original in list(output['muscles']):
        mirrored=copy.deepcopy(original); mirrored['name']=original['name'][:-1]+'l'
        mirrored['mirrored_from']=original['name']
        mirrored['source_neutral_moment_arms_m']={k.replace('_r','_l'):v for k,v in original['source_neutral_moment_arms_m'].items() if not k.endswith('_l')}
        for point in mirrored['points']:
            for field in ('name','parent_frame','body','condition_coordinate'):
                if point.get(field): point[field]=point[field].replace('_r','_l')
            for field in ('location_parent_m','location_ground_m'):
                point[field][2]*=-1
            point['source_xml']=point['source_xml'].replace('_r','_l')
            for axis,function in enumerate(point.get('motion_functions',[])):
                function['coordinate']=function['coordinate'].replace('_r','_l')
                if axis==2:
                    function['coefficients']=[[-x for x in c] for c in function['coefficients']]
                    function['endpoints']=[[-x for x in c] for c in function['endpoints']]
        output['muscles'].append(mirrored)
    mapping=json.loads((REPO/'Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json').read_text())
    output['actuator_map_source']='Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json'
    output['actuators_27']=mapping['actuators']
    (HERE/'robotbody_source_inventory.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
    types={}
    for muscle in output['muscles']:
        for point in muscle['points']:
            types[point['type']]=types.get(point['type'],0)+1
    print(json.dumps(dict(muscles=len(output['muscles']),point_types=types,
                         muscles_with_wraps=sum(x['path_wrap_count']>0 for x in output['muscles']),
                         joint_names=list(output['joints'])),indent=2))


if __name__ == '__main__':
    main()
