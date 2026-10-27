"""Weld the supplied resolved CAD foot geometry to each distal shank."""
import json
import xml.etree.ElementTree as ET
import numpy as np
from add_torso import values

def distal_center(side,cad,here):
    component=next(x for x in cad['components'] if ('05_01_TI_R' if side=='r' else '05_02_TI') in x['name'])
    t=np.array(component['transform']); r=t[:9].reshape(3,3).T
    raw=(here/'ah001_cad_meshes/meshes'/('shank_'+side+'.STL')).read_bytes()
    n=int(np.frombuffer(raw,'<u4',count=1,offset=80)[0])
    triangles=np.frombuffer(raw,dtype=np.dtype([('normal','<f4',(3,)),('vertices','<f4',(3,3)),('attr','<u2')]),count=n,offset=84)
    local=(triangles['vertices'].reshape(-1,3)-t[9:12])@r
    end=local[local[:,1]<local[:,1].min()+.0005]
    c=(end.min(axis=0)+end.max(axis=0))/2; c[1]=local[:,1].min()
    return r@c+t[9:12]

def add_feet(root,bodies,origins,properties,audit,here):
    from build_ah001 import inertia_for_group,localize_mesh,cad_to_isb
    foot=json.loads((audit/'foot_inventory.json').read_text())
    cad=json.loads((audit/'robot_mass_geometry.json').read_text())
    components={c['name']:c for c in foot['components']}
    assert len(components)==22 and all(c['mass_properties_density_one'] for c in components.values())
    densities={n:7850. if 'axle' in n.lower() and 'lock' not in n.lower() else 1200. for n in components}
    mass,center,tensor=inertia_for_group(list(components),components,densities,np.eye(3),np.zeros(3))
    raw=(audit/'foot_cad.stl').read_bytes(); n=int(np.frombuffer(raw,'<u4',count=1,offset=80)[0])
    vertices=np.frombuffer(raw,dtype=np.dtype([('normal','<f4',(3,)),('vertices','<f4',(3,3)),('attr','<u2')]),count=n,offset=84)['vertices'].reshape(-1,3)
    exterior=components['Ankle Exterior-1']; t=np.array(exterior['transform']);r=t[:9].reshape(3,3).T
    extcenter=r@np.array(exterior['mass_properties_density_one'][0][:3])+t[9:12]
    mount=np.array([extcenter[0],vertices[:,1].max(),extcenter[2]])
    result={}
    for side in ('r','l'):
        name='foot_'+side
        rotation=np.array([[0.,0,1],[0,1,0],[-1,0,0]])
        if side=='l': rotation=np.diag([1,1,-1])@rotation
        anchor=cad_to_isb(distal_center(side,cad,here))
        origins[name]=anchor
        body=ET.SubElement(bodies['shank_'+side],'body',name=name,pos=values(anchor-origins['shank_'+side]))
        bodies[name]=body
        com=rotation@(center-mount); inertia=rotation@tensor@rotation.T
        ET.SubElement(body,'inertial',mass=str(mass),pos=values(com),
          fullinertia=values([inertia[0,0],inertia[1,1],inertia[2,2],inertia[0,1],inertia[0,2],inertia[1,2]]))
        bounds,files=localize_mesh(audit/'foot_cad.stl',here/'meshes'/(name+'.stl'),rotation@mount,rotation)
        world_bounds=(np.array(bounds)-rotation@mount+anchor).tolist()
        for index,filename in enumerate(files):
            mesh=name+'_%02d'%index
            ET.SubElement(root.find('asset'),'mesh',name=mesh,file=filename)
            ET.SubElement(body,'geom',name=mesh+'_cad',type='mesh',mesh=mesh,contype='1',conaffinity='2',friction='.9 .005 .0001',rgba='.3 .35 .4 1',condim='3')
        prop=dict(mass_kg=mass,center_of_mass=(anchor+com).tolist(),inertia_com=inertia.tolist(),mesh_bounds=world_bounds)
        properties[name]=prop
        result[side]=dict(**prop,mount_source_m=mount.tolist(),mount_pelvis_frame_m=anchor.tolist(),
                         source=foot['source'],attachment='rigid weld at distal shank; all foot internal joints welded',
                         density_assumption='solid Onyx; axle shafts steel, collars and locks Onyx',
                         mirrored=side=='l')
    return result
