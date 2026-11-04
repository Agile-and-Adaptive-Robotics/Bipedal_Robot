"""Catalog-envelope asset library; proposed brackets are explicitly labelled."""
from pathlib import Path
import bpy,json,math
from mathutils import Vector
P=Path(__file__).resolve().parents[1]/'packaging'
D=json.loads((P/'hardware/hardware_candidates.json').read_text())
bpy.ops.object.select_all(action='SELECT');bpy.ops.object.delete(use_global=False)
s=bpy.context.scene;s.unit_settings.system='METRIC';s.world.color=(1,1,1)
s.render.engine='CYCLES';s.cycles.samples=16;s.render.threads_mode='FIXED';s.render.threads=4
s.render.resolution_x=1600;s.render.resolution_y=1000;s.render.resolution_percentage=100
s.view_settings.view_transform='Standard'
try:s.view_settings.look='None'
except TypeError:pass
def mat(name,color,metal=0):
 m=bpy.data.materials.new(name);m.diffuse_color=(*color,1);m.use_nodes=True
 bs=m.node_tree.nodes.get('Principled BSDF');bs.inputs['Base Color'].default_value=(*color,1);bs.inputs['Metallic'].default_value=metal;bs.inputs['Roughness'].default_value=.45;return m
wheel=mat('Catalog sheave envelope',(.08,.28,.40),.4)
proposal=mat('Proposed bracket - not supplier geometry',(.75,.34,.09),.3)
steel=mat('Pin envelope - structural sizing pending',(.3,.34,.36),.7)
black=mat('Housing',(.035,.035,.04));textmat=mat('Labels',(.03,.03,.04))
current=None
def own(o,name,material):
 o.name=name
 for c in list(o.users_collection):c.objects.unlink(o)
 current.objects.link(o);o.data.materials.append(material);return o
def cyl(name,r,depth,pos,material):
 bpy.ops.mesh.primitive_cylinder_add(vertices=64,radius=r,depth=depth,location=pos);return own(bpy.context.object,name,material)
def box(name,dims,pos):
 bpy.ops.mesh.primitive_cube_add(size=1,location=pos);o=own(bpy.context.object,name,proposal);o.dimensions=dims;return o
def label(body,x,y,size=.005):
 curve=bpy.data.curves.new('Label','FONT');curve.body=body;curve.size=size
 o=bpy.data.objects.new(body,curve);current.objects.link(o);o.location=(x,y,.025);o.data.materials.append(textmat)
def sheave(f,pos):
 # Revolve a closed bore-to-rim profile. Groove is schematic; catalog OD,
 # width and bore are exact, but the catalog pulley diameter is not a
 # verified wire pitch diameter or a manufacturing groove drawing.
 r=f['od_mm']/2000;w=f['width_mm']/1000;b=f['bore_mm']/2000
 rg=f['catalog_pulley_diameter_mm']/2000
 profile=[(b,-w/2),(r,-w/2),(r,-w*.30),(rg,0),(r,w*.30),(r,w/2),(b,w/2)]
 vertices=[];faces=[];N=96
 for rad,z in profile:
  vertices.extend((pos[0]+rad*math.cos(i*2*math.pi/N),pos[1]+rad*math.sin(i*2*math.pi/N),pos[2]+z) for i in range(N))
 for j in range(len(profile)):
  for i in range(N):faces.append((j*N+i,j*N+(i+1)%N,((j+1)%len(profile))*N+(i+1)%N,((j+1)%len(profile))*N+i))
 me=bpy.data.meshes.new('Catalog envelope mesh');me.from_pydata(vertices,[],faces);me.update()
 o=bpy.data.objects.new(f['sheave_part']+' sheave envelope',me);current.objects.link(o);o.data.materials.append(wheel)
 for k,v in f.items():o[k]=str(v)
 o['geometry_status']='Verified catalog OD/width/bore. Schematic groove; not supplier CAD.'
 return o
for index,f in enumerate(D['families']):
 current=bpy.data.collections.new(f['id']+' candidate family - capacity of assembly unverified');s.collection.children.link(current)
 x=(index-1)*.16;diam=f['od_mm']/1000;width=f['width_mm']/1000
 # A loose sheave and a two-sheave block envelope use the same part family.
 sheave(f,(x,.075,0))
 label(f"{f['id']}: McMaster {f['sheave_part']}",x-.065,.123)
 label(f"OD {f['od_mm']:g} / width {f['width_mm']:g} mm",x-.065,.113,.004)
 label(f"Catalog capacity {f['catalog_capacity_lbf']} lbf",x-.065,.103,.004)
 plate=f['proposed_plate_thickness_mm']/1000;gap=.001
 zcenters=[-(width+gap)/2,(width+gap)/2]
 for z in zcenters:sheave(f,(x,0,z))
 stack=2*width+gap
 for sign in (-1,1):
  o=box('PROPOSED side plate', (diam+.004,diam+.004,plate),(x,0,sign*(stack/2+gap+plate/2)))
  o['status']='Illustrative packaging allowance, not a rated or purchased bracket.'
 cyl('Unrated axle envelope',f['bore_mm']/2000,stack+2*gap+2*plate+.008,(x,0,0),steel)
 label('Two-sheave block ENVELOPE',x-.065,-.057,.004)
 label('Reeving and mount design pending',x-.065,-.067,.0035)
current=bpy.data.collections.new('Bowden dimensional references - loads unverified');s.collection.children.link(current)
for x,od,caption in [(-.13,.005,'Bicycle: 5 mm housing / 6 mm ferrule'),(.06,.007747,'Sava C277PL: 7.747 mm finished OD')]:
 cyl(caption,od/2,.12,(x,-.145,0),black)
 if od==.005:
  for sign in (-1,1):
   o=cyl('BOT192-BJ envelope: length assumed 12 mm',.003,.012,(x,-.145,sign*.061),steel)
   o['catalog_verified']='6 mm OD and 5 mm housing fit';o['assumed_length_mm']=12
 label(caption,x-.06,-.18,.0036)
current=bpy.data.collections.new('Scope');s.collection.children.link(current)
label('Catalog envelopes + proposed block allowances. NOT released load-rated assemblies.',-.23,-.215,.004)
label('Three hardware families; reeving, cable rating, stops, and cyclic life still require verification.',-.23,-.225,.0037)
s['catalog_and_assumptions']=json.dumps(D)
bpy.data.texts.new('CATALOG SOURCES AND ASSUMPTIONS').write(json.dumps(D,indent=2))
camdata=bpy.data.cameras.new('Camera');cam=bpy.data.objects.new('Camera',camdata);s.collection.objects.link(cam)
cam.location=(0,-.32,1.0);target=Vector((0,-.04,0));cam.rotation_euler=(target-cam.location).to_track_quat('-Z','Y').to_euler();camdata.type='ORTHO';camdata.ortho_scale=.76;s.camera=cam
ld=bpy.data.lights.new('Area','AREA');ld.energy=3;ld.size=.6;lo=bpy.data.objects.new('Area',ld);s.collection.objects.link(lo);lo.location=(0,0,.7)
s.world.use_nodes=True;s.world.node_tree.nodes['Background'].inputs['Color'].default_value=(1,1,1,1)
s.world.node_tree.nodes['Background'].inputs['Strength'].default_value=.8
for a in bpy.context.screen.areas:
 if a.type=='VIEW_3D':a.spaces.active.region_3d.view_perspective='CAMERA'
path=P/'output/AH001_Hardware_Candidates.blend'
bpy.ops.wm.save_as_mainfile(filepath=str(path),compress=True)
s.render.filepath=str(P/'output/AH001_Hardware_Candidates.png');bpy.ops.render.render(write_still=True)
print('Saved',path)
