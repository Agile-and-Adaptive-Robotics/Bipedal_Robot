"""Offline CAD/KompareXXVII inspection scene. Reads saved designs, never solves them.

The CAD uses the audited SOLIDWORKS pelvis datum and rigid groups, not a fitted
registration. Survey-to-CAD attachment registration remains a preview assumption.
"""
from pathlib import Path
import base64, csv, gzip, hashlib, json, math, struct, re
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation, Slerp
from scipy.spatial import cKDTree
import trimesh

HERE=Path(__file__).resolve().parent
REP=HERE.parent
REPO=REP.parents[1]
OUT=REP/'packaging/output'; OUT.mkdir(parents=True,exist_ok=True)
RT=REPO/'Code/MuJoCo_SNS/robot_torso'
AUD=RT.parent/'spinal/reports_20261027_robot_control'
PKG=REPO/'Code/Matlab/Mesh_Optimization/Results/KompareXXVII'
NS={'s':'urn:mathworks:SimscapeMultibody:import'}
partdata={}; registrations={}; stl_cache={}
for folder,filename in [('export','ah001_full.xml.xml'),('export_foot','foot_ankle.xml.xml')]:
    xr=ET.parse(REP/'cad'/folder/filename).getroot()
    for part in xr.findall('.//s:Part',NS):
        cm=part.find('s:MassProperties/s:CenterOfMass',NS); geo=part.find('s:GeometryFile',NS)
        if cm is not None and geo is not None:
            partdata[part.get('name')]=(np.fromstring(cm.text,sep=' '),REP/'cad'/folder/geo.get('name'))

def pose(R=np.eye(3),t=np.zeros(3)):
    a=np.eye(4); a[:3,:3]=R; a[:3,3]=t; return a
def apply(T,P): return np.asarray(P)@T[:3,:3].T+T[:3,3]
def stl(path):
    if path in stl_cache:return stl_cache[path]
    raw=path.read_bytes(); n=struct.unpack_from('<I',raw,80)[0]
    assert len(raw)==84+50*n,path
    d=np.frombuffer(raw,dtype=np.dtype([('n','<f4',3),('v','<f4',(3,3)),('a','<u2')]),offset=84,count=n)
    # SOLIDWORKS SaveAs STL coordinates are millimetres in these files.
    tri=d['v'].astype(float)/1000.0
    # These STL exports also moved each part into the positive octant. Recover
    # the lost translation by matching actual STEP points, seeded by volume COM.
    if path.stem not in partdata: raise ValueError('Missing local CAD datum: '+str(path))
    cm,step=partdata[path.stem]
    m=trimesh.Trimesh(vertices=tri.reshape(-1,3),faces=np.arange(tri.size//3).reshape(-1,3),process=True)
    seed=cm-m.center_mass
    steptext=step.read_text(errors='replace')
    scale=.0254 if re.search(r"CONVERSION_BASED_UNIT\s*\(\s*'INCH'",steptext) else .001
    pts=np.array([list(map(float,s.split(','))) for s in re.findall(r'CARTESIAN_POINT\s*\(\s*[^,]*,\s*\(([^)]+)\)',steptext,re.S)])*scale
    v=np.unique(tri.reshape(-1,3),axis=0); dist,ix=cKDTree(pts).query(v+seed)
    good=dist<.00015
    if good.sum()<4:raise ValueError('Cannot recover STL offset: '+str(path))
    offsets=pts[ix[good]]-v[good]
    rounded=np.round(offsets,6); unique,counts=np.unique(rounded,axis=0,return_counts=True)
    cluster=offsets[np.all(rounded==unique[counts.argmax()],axis=1)]
    shift=np.median(cluster,axis=0)
    exact_dist=cKDTree(pts).query(v+shift)[0]
    matches=int((exact_dist<2e-7).sum())
    if matches<4:raise ValueError('STEP vertex registration failed: '+str(path))
    registrations[path.stem]=dict(offset_m=shift.tolist(),matching_STEP_vertices=matches,
        max_cluster_spread_m=float(np.max(np.linalg.norm(cluster-shift,axis=1))),seed_COM_difference_m=float(np.linalg.norm(shift-seed)))
    tri+=shift;stl_cache[path]=tri
    return tri
def clean(x):
    if isinstance(x,np.ndarray): return clean(x.tolist())
    if isinstance(x,dict): return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)) and not np.isfinite(x): return None
    if isinstance(x,np.generic): return x.item()
    return x

root=ET.parse(RT/'ah001_neural.xml').getroot(); tree={}
def walk(el,parent=None):
    for b in el.findall('body'):
        name=b.get('name'); q=np.fromstring(b.get('quat','1 0 0 0'),sep=' ')
        T=pose(Rotation.from_quat(np.r_[q[1:],q[0]]).as_matrix(),np.fromstring(b.get('pos','0 0 0'),sep=' '))
        j=b.find('joint'); axis=np.fromstring(j.get('axis','0 0 1'),sep=' ') if j is not None else np.array([0.,0.,1.])
        tree[name]=dict(parent=parent,T=T,axis=axis)
        walk(b,name)
walk(root.find('worldbody'))
def fk(q):
    W={}
    for b,v in tree.items():
        T=v['T'].copy()
        if b in q: T=T@pose(Rotation.from_rotvec(v['axis']*q[b]).as_matrix())
        W[b]=(W[v['parent']] if v['parent'] else np.eye(4))@T
    return W
W0=fk({}); ph=W0['pelvis'][1,3]
closures=json.loads((RT/'ah001_loop_closures.json').read_text())
anchors={}
for c in closures:
    b1,b2=c['body1'],c['body2']; a=np.array(c['anchor_body1_m'])
    anchors[b1[-1]]=(a,apply(np.linalg.inv(W0[b2]),apply(W0[b1],a)))

# Exact CAD loop continuation, prescribed shank orientation change.
cadposes=[]; loop_res=[]; prev={s:np.zeros(3) for s in ('r','l')}
for deg in np.linspace(0,-90,19):
    q={}
    for s in ('r','l'):
        a,b=anchors[s]
        def residual(x):
            w=fk({f'back_{s}':x[0],f'shank_{s}':x[1],f'front_{s}':x[2]})
            dp=apply(w[f'front_{s}'],a)-apply(w[f'shank_{s}'],b)
            return np.r_[dp,x[0]+x[1]+math.radians(deg)]
        sol=least_squares(residual,prev[s],xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=300)
        prev[s]=sol.x
        err=float(np.linalg.norm(residual(sol.x)[:3])); loop_res.append(err)
        if err>1e-5: raise RuntimeError(f'Closed-loop CAD pose failed at {deg}: {err}')
        q.update({f'back_{s}':sol.x[0],f'shank_{s}':sol.x[1],f'front_{s}':sol.x[2]})
    cadposes.append(dict(knee_deg=float(deg),W=clean(fk(q))))

# SOLIDWORKS component local geometry -> audited anatomical pelvis -> body local.
cad=json.loads((AUD/'robot_mass_geometry.json').read_text())
audit=json.loads((AUD/'robot_geometry_audit.json').read_text())
pc=next(c for c in cad['components'] if '/02_01_PE_' in c['name'])
pt=np.array(pc['transform']); CR=pt[:9].reshape(3,3); CT=pt[9:12]
group={n:g for g,ns in audit['rigid_groups'].items() for n in ns}
meshes=[]; glb=trimesh.Scene()
def addmesh(name,body,tri):
    # Deduplicate without dropping triangles or opening the CAD surfaces.
    m=trimesh.Trimesh(vertices=tri.reshape(-1,3),faces=np.arange(tri.size//3).reshape(-1,3),process=True)
    v=m.vertices.astype('<f4'); f=m.faces.astype('<u4')
    meshes.append(dict(name=name,body=body,v=base64.b64encode(v.tobytes()).decode(),f=base64.b64encode(f.tobytes()).decode()))
    cp=m.copy(); cp.apply_transform(W0[body]); cp.visual.face_colors=[120,135,146,255]
    glb.add_geometry(cp,node_name=name,geom_name=name)
for c in cad['components']:
    name=c['name']; stem=name.rsplit('/',1)[-1].rsplit('-',1)[0]
    file=REP/'cad/stl'/f'{stem}.STL'
    if name not in group or not file.exists(): continue
    t=np.array(c['transform']); R=CR@t[:9].reshape(3,3).T; p=CR@(t[9:12]-CT)
    P=stl(file)@R.T+p+np.array([0,ph,0]); body=group[name]
    P=apply(np.linalg.inv(W0[body]),P)
    addmesh(name,body,P)

# Feet use the exported assembly's exact relative part placements. The foot-to-
# shank mount is inherited from ah001_neural, not a new measured mount.
ns={'s':'urn:mathworks:SimscapeMultibody:import'}
fx=ET.parse(REP/'cad/export_foot/foot_ankle.xml.xml').getroot()
fp={}
for el in fx.findall('.//s:Instance',ns):
    tr=el.find('s:Transform',ns)
    if tr is None: continue
    fp[el.get('name')]=pose(np.fromstring(tr.find('s:Rotation',ns).text,sep=' ').reshape(3,3),np.fromstring(tr.find('s:Translation',ns).text,sep=' '))
foot_inventory=json.loads((AUD/'foot_inventory.json').read_text())
foot_mounts=json.loads((RT/'ah001_build_audit.json').read_text())['feet']
foot_instance={}
for c in foot_inventory['components']:
    t=np.array(c['transform']);foot_instance[c['name']]=pose(t[:9].reshape(3,3).T,t[9:12])
for name,T in foot_instance.items():
    stem=name.rsplit('-',1)[0]; file=REP/'cad/stl'/f'{stem}.STL'
    if not file.exists(): continue
    p=apply(T,stl(file))
    for s in ('r','l'):
        # The assembly has width/up/anterior axes. The MuJoCo foot body uses
        # anatomical anterior/up/right axes and its origin is the distal mount,
        # not the Ankle Exterior part datum. Match foot_geometry.py verbatim.
        Rfoot=np.array([[0.,0,1],[0,1,0],[-1,0,0]])
        if s=='l':Rfoot=np.diag([1,1,-1])@Rfoot
        pp=(p-np.array(foot_mounts[s]['mount_source_m']))@Rfoot.T
        if np.linalg.det(Rfoot)<0:pp=pp[:,[0,2,1],:]
        addmesh('foot_'+s+'/'+name,'foot_'+s,pp)

spine=trimesh.creation.box(extents=[.02,.25,.02]);spine.apply_translation([0,.125,0])
addmesh('Rigid 0.25 m spine','lumbar_torso',spine.triangles)

routes=json.loads((REP/'packaging/build/saved_routes.json').read_text())
configs={r['site']:r for r in csv.DictReader((PKG/'route_configurations.csv').open()) if r['side']=='selected'}
def body(fr):
    if fr.startswith('pelvis'):return 'pelvis'
    if fr.startswith('femur'):return 'thigh_r'
    if fr.startswith('tibia'):return 'shank_r'
    if fr.startswith('torso'):return 'lumbar_torso'
    return 'foot_r'
def interpT(ts,angles,x):
    if np.ptp(angles)<1e-10: return ts[:,:,0]
    x=float(np.clip(x,angles[0],angles[-1])); R=Slerp(angles,Rotation.from_matrix(np.moveaxis(ts[:3,:3,:],2,0)))([x]).as_matrix()[0]
    return pose(R,np.array([np.interp(x,angles,ts[k,3,:]) for k in range(3)]))
result=[]
for r in routes:
    sp=r['spec']; frames=sp['frames']; P=np.array(sp['points']); cross=np.array(sp['cross']); T=np.array(r['T'])
    ranges=np.array(sp['kin']['thetaRangesDeg']); n1,n2=np.array(r['sContraction']).shape
    a1=np.linspace(*ranges[0],n1); a2=np.linspace(*ranges[1],n2)
    Ts=[interpT(T[:,:,:,0],a1,0),interpT(T[:,:,:,1],a2,0)]
    chain=[np.eye(4),Ts[0],Ts[0]@Ts[1]]
    pointBodies=[]; local=[]
    for i,p in enumerate(P):
        f=0 if i+1<cross[0] else (2 if cross[0]==cross[1] or i+1>=cross[1] else 1)
        b=body(frames[f]); b0=body(frames[0]); neutral=apply(W0[b0]@chain[f],p)
        pointBodies.append(b); local.append(apply(np.linalg.inv(W0[b]),neutral))
    cfg=configs[r['site']]
    offsets=np.array(sp['bundleOffsets'],dtype=float).reshape(-1,3)
    if len(offsets)==0: offsets=np.zeros((1,3))
    item=dict(site=r['site'],name=cfg['muscle'],bodies=pointBodies,local=local,offsets=offsets,offset_body=body(frames[0]),
        rest=r['RestingL'],diameter=r['Diameter']/1000,minimum=r['Kmax'],fit=r['FittingLength'],tendon=r['TendonL'],
        n=int(cfg['nBPA']),gain=float(cfg['G_survey']),ranges=ranges,angles=[a1,a2],dofs=sp['kin']['dofs'],
        contraction=r['sContraction'],reported_strain=r['saved_strain'],slack=r['PulleySlack'],infeasible=r['PulleyInfeasible'],source=r['source'],
        note=sp['notes'],frames=frames,cross=cross,neutral_chain=chain)
    result.append(item)

data=dict(meshes=meshes,poses=cadposes,routes=result,W0=W0,frame='x anterior, y superior, z right',
    warning='Attachment registration is a datum-based preview. Physical tackle/sheave layouts are not saved. Overlaps are screening candidates, not certified clearances.',
    audit=dict(route_count=len(result),cad_mesh_count=len(meshes),loop_max_error_m=max(loop_res),
        registration='Audited CAD pelvis normalization; source survey chains folded exactly at zero degrees, then attached to preview CAD bodies. No fitted global CAD transform.',
        diameter='Calibrated inextensible braid endpoints; saved sContraction sets actuator length. Values outside calibration are flagged.',
        output='Saved-route envelopes, not a newly optimized or validated physical installation.'))
blob=gzip.compress(json.dumps(clean(data),separators=(',',':'),allow_nan=False).encode())
template=(HERE/'packaging_view.html').read_text(encoding='utf-8')
(OUT/'AH001_Packaging.html').write_text(template.replace('__DATA__',base64.b64encode(blob).decode()),encoding='utf-8')
(OUT/'scene_data.json').write_text(json.dumps(clean(data),separators=(',',':')),encoding='utf-8')
(OUT/'verification.json').write_text(json.dumps(clean(data['audit']),indent=2))
(OUT/'mesh_registration.json').write_text(json.dumps(registrations,indent=2))
print(json.dumps(data['audit'],indent=2))
# CAD-only GLB is a faithful editable import. BPA geometry is exported by viewer.
glb.export(OUT/'AH001_CAD.blend_import.glb')
print('HTML and GLB:',OUT)
