"""Fixed-length schematic cable driving the existing crossed knee linkage.

Coordinates are from the source figure, in schematic units. This is a
kinematic illustration, not a dynamic result or a new robot design.
"""
import numpy as np
import math
A=np.array([6.25,3.85]); B=np.array([6.65,4.05])
C0=np.array([6.75,3.2]); D0=np.array([7.05,3.35])
P0=np.array([7.8,2.35]); O=np.array([7.35,5.1]); R=.45
def rot(q):return np.array([[math.cos(q),-math.sin(q)],[math.sin(q),math.cos(q)]])
def pose(q):
 D=A+rot(q)@(D0-A)
 # C lies on circles about fixed B and moving D. Retain the source branch.
 r=np.linalg.norm(C0-B); s=np.linalg.norm(C0-D0); v=D-B; d=np.linalg.norm(v)
 a=(r*r-s*s+d*d)/(2*d); h2=r*r-a*a
 if h2 < -1e-10:raise ValueError('Four-bar closure unavailable')
 u=v/d; perp=np.array([-u[1],u[0]]); base=B+a*u
 cs=[base+np.sqrt(max(0,h2))*perp,base-np.sqrt(max(0,h2))*perp]
 # The lower intersection continues the initial crossed-linkage branch.
 C=min(cs,key=lambda x:x[1])
 theta=math.atan2(*(C-D)[::-1])-math.atan2(*(C0-D0)[::-1])
 Q=rot(theta); t=D-Q@D0; P=Q@P0+t
 vec=P-O; dist=np.linalg.norm(vec)
 if dist<=R:raise ValueError('Output attachment inside pulley')
 tangent_angle=math.atan2(vec[1],vec[0])+math.acos(R/dist)
 T=O+R*np.array([math.cos(tangent_angle),math.sin(tangent_angle)])
 out=R*(math.pi/2-tangent_angle)+np.linalg.norm(P-T)
 lam=np.linalg.solve(np.column_stack([D-A,-(C-B)]),B-A)[0]
 I=A+lam*(D-A)
 return dict(D=D,C=C,P=P,I=I,Q=Q,t=t,theta=theta,tangent=T,angle=tangent_angle,out=out)
L0=pose(0)['out']
def contraction(q):return (L0-pose(q)['out'])/2
def end_angle(target=.42):
 prev=0.
 for q in np.linspace(.001,.65,300):
  try:d=contraction(q)
  except ValueError:break
  if d>=target:
   lo,hi=prev,q
   for _ in range(45):
    mid=(lo+hi)/2
    if contraction(mid)<target:lo=mid
    else:hi=mid
   return (lo+hi)/2
  prev=q
 raise ValueError('Requested shortening exceeds the source linkage range')
if __name__=='__main__':
 qmax=end_angle();errors=[]
 for q in np.linspace(0,qmax,101):
  p=pose(q);d=contraction(q)
  errors.append(abs((2*d+p['out'])-L0))
  assert abs(np.linalg.norm(p['D']-A)-np.linalg.norm(D0-A))<1e-10
  assert abs(np.linalg.norm(p['C']-B)-np.linalg.norm(C0-B))<1e-10
  assert abs(np.linalg.norm(p['C']-p['D'])-np.linalg.norm(C0-D0))<1e-10
 print('Drive angle deg:',np.degrees(qmax),'BPA travel:',contraction(qmax),'output cable travel:',L0-pose(qmax)['out'],'max cable length residual:',max(errors))
