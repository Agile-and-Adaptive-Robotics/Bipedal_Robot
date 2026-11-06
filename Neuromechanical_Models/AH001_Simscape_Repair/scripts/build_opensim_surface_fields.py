"""Signed surface fields from the exact displayed OpenSim geometry."""
from pathlib import Path
import json,base64
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy,numpy_to_vtk
from scipy.io import savemat
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/audit_20261105'
D=json.loads((OUT/'opensim_verified_scene.json').read_text());dest=OUT/'surface_fields';dest.mkdir(exist_ok=True)
reports=[]
for body in ('pelvis','femur_r','tibia_r','talus_r','calcn_r','toes_r','torso'):
 meshes=[m for m in D['meshes'] if m['body']==body];poly=[];vertices=[];open_edges=0
 for m in meshes:
  v=np.frombuffer(base64.b64decode(m['vertices']),dtype='<f4').reshape(-1,3).astype(float)
  f=np.frombuffer(base64.b64decode(m['faces']),dtype='<u4').reshape(-1,3).astype(np.int64);vertices.append(v)
  pts=vtk.vtkPoints();pts.SetData(numpy_to_vtk(v,deep=True));cells=vtk.vtkCellArray();cells.SetCells(len(f),numpy_to_vtk(np.column_stack([np.full(len(f),3),f]).ravel(),deep=True,array_type=vtk.VTK_ID_TYPE))
  p=vtk.vtkPolyData();p.SetPoints(pts);p.SetPolys(cells)
  clean=vtk.vtkCleanPolyData();clean.SetInputData(p);clean.Update()
  edges=vtk.vtkFeatureEdges();edges.SetInputData(clean.GetOutput());edges.BoundaryEdgesOn();edges.NonManifoldEdgesOn();edges.FeatureEdgesOff();edges.ManifoldEdgesOff();edges.Update();open_edges+=edges.GetOutput().GetNumberOfCells()
  normals=vtk.vtkPolyDataNormals();normals.SetInputData(clean.GetOutput());normals.ConsistencyOn();normals.AutoOrientNormalsOn();normals.SplittingOff();normals.Update();poly.append(normals.GetOutput())
 if open_edges:
  reports.append(dict(body=body,status='OPEN SURFACE - NOT USED',boundary_or_nonmanifold_edges=open_edges));print(reports[-1],flush=True);continue
 V=np.vstack(vertices);lo=V.min(0)-.05;hi=V.max(0)+.05;step=.003
 count=np.ceil((hi-lo)/step).astype(int)+1;axes=[np.linspace(lo[k],hi[k],count[k]) for k in range(3)];sdf=None
 for p in poly:
  distance=vtk.vtkImplicitPolyDataDistance();distance.SetInput(p)
  sample=vtk.vtkSampleFunction();sample.SetImplicitFunction(distance);sample.SetModelBounds(*[a for k in range(3) for a in (lo[k],hi[k])]);sample.SetSampleDimensions(*map(int,count));sample.ComputeNormalsOff();sample.Update()
  values=vtk_to_numpy(sample.GetOutput().GetPointData().GetScalars()).reshape(count[2],count[1],count[0]).transpose(2,1,0)
  sdf=values if sdf is None else np.minimum(sdf,values)
 uncertainty=float(np.linalg.norm([a[1]-a[0] for a in axes]))
 savemat(dest/(body+'.mat'),dict(x=axes[0],y=axes[1],z=axes[2],distance=sdf.astype(np.float32),uncertainty_m=uncertainty,mesh_bbox=np.vstack([V.min(0),V.max(0)]),body=body),do_compression=True)
 reports.append(dict(body=body,status='CLOSED SURFACE',triangles=sum(p.GetNumberOfCells() for p in poly),grid=count.tolist(),uncertainty_m=uncertainty,negative_voxels=int(np.sum(sdf<0)),boundary_or_nonmanifold_edges=open_edges))
 print(reports[-1],flush=True)
(dest/'manifest.json').write_text(json.dumps(dict(source=D['source'],convention='Negative inside, per-piece union; query subtracts full cell diagonal as a conservative interpolation allowance. These fields describe OpenSim geometry, not a guessed CAD registration.',bodies=reports),indent=2))
