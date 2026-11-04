"""Shared metric preview geometry for corrected GLB and native Blender scenes."""
import math
import numpy as np

def apply(T,p):
 T=np.asarray(T);return np.asarray(p)@T[:3,:3].T+T[:3,3]
def unit(v):return np.asarray(v)/max(np.linalg.norm(v),1e-12)
def contraction(r,deg):
 k=next((i for i,v in enumerate(r['dofs']) if 'knee' in v),None)
 if k is None:
  i=int(np.abs(r['angles'][0]).argmin());j=int(np.abs(r['angles'][1]).argmin());return r['contraction'][i][j]
 other=int(np.abs(r['angles'][1-k]).argmin());s=np.asarray(r['contraction']);s=s[:,other] if k==0 else s[other,:]
 return float(np.interp(deg,r['angles'][k],s))
def diameter(r,L):
 # All replacement actuators are 20 mm. This law applies independently to
 # each half of a split BPA, which preserves the same contraction fraction.
 rest=r['rest'];lm=.75*rest;pn=math.sqrt((rest**2-lm**2)/(.038**2-.020**2))
 L=np.clip(L,lm,rest)
 return math.sqrt(rest**2+(pn*.020)**2-L**2)/pn
def sample_path(P,start,end,M=32):
 P=np.asarray(P);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(P,axis=0),axis=1))]
 keep=np.r_[True,np.diff(arc)>1e-10];arc=arc[keep];P=P[keep]
 if len(P)<2:return np.repeat(P[:1],M,axis=0)
 ds=np.linspace(np.clip(start,0,arc[-1]),np.clip(end,0,arc[-1]),M)
 return np.stack([np.interp(ds,arc,P[:,k]) for k in range(3)],axis=1)
def between(P,a,b):
 arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(P,axis=0),axis=1))]
 return P[(arc>a)&(arc<b)]
def tube(P,r,N=16):
 P=np.asarray(P);v=[]
 for i,p in enumerate(P):
  d=unit(P[min(i+1,len(P)-1)]-P[max(0,i-1)])
  ref=np.array([0,1,0]) if abs(d[1])<.9 else np.array([1,0,0]);u=unit(np.cross(d,ref));w=np.cross(d,u)
  v.extend(p+r*(math.cos(j*2*math.pi/N)*u+math.sin(j*2*math.pi/N)*w) for j in range(N))
 f=[]
 for i in range(1,len(P)):
  for j in range(N):f.append(((i-1)*N+j,(i-1)*N+(j+1)%N,i*N+(j+1)%N,i*N+j))
 f.extend([tuple(range(N-1,-1,-1)),tuple((len(P)-1)*N+j for j in range(N))])
 return np.array(v),f
def parts(r,W,deg):
 """Return fixed-topology components and explicit shared insertion pins."""
 P=np.array([apply(W[b],p) for b,p in zip(r['bodies'],r['local'])]);total=np.linalg.norm(np.diff(P,axis=0),axis=1).sum()
 L=r['rest']-contraction(r,deg);D=diameter(r,L);fit=r['fit'];out={}
 def add(name,line,rad,kind,M=32):
  line=sample_path(np.asarray(line),0,np.linalg.norm(np.diff(line,axis=0),axis=1).sum(),M)
  v,f=tube(line,rad);out[name]=dict(vertices=v,faces=f,kind=kind)
 # Saved wrap-elimination states can repeat terminal points. Find a genuine
 # segment before orienting a pinned sleeve; a zero tangent collapses it.
 forward=next((p-P[0] for p in P[1:] if np.linalg.norm(p-P[0])>1e-9),None)
 backward=next((p-P[-1] for p in P[-2::-1] if np.linalg.norm(p-P[-1])>1e-9),None)
 if forward is None or backward is None:raise ValueError(r['site']+': coincident route points')
 u=unit(forward);v=unit(backward)
 if r['split_bowden']:
  A=P[0]+u*(L/2+2*fit);B=P[-1]+v*(L/2+2*fit)
  middle=between(P,L/2+2*fit,total-(L/2+2*fit))
  housing=np.vstack([A,middle,B])
  for bi,offset in enumerate(r['offsets'],1):
   off=np.asarray(W[r['offset_body']])[:3,:3]@offset
   for label,O,d in [('proximal',P[0],u),('distal',P[-1],v)]:
    # Each end bundle follows its own pinned sleeve orientation. An origin
    # frame offset must not slide distal actuators along their length.
    transverse=off-d*np.dot(off,d)
    x=O+transverse+d*fit;y=x+d*(L/2);end=y+d*fit
    add(f'BPA {bi}/{label} sleeve',[x,y],D/2,'sleeve',24)
    add(f'BPA {bi}/{label} outer fitting',[O,x],.0106,'fitting',3)
    add(f'BPA {bi}/{label} inner fitting',[y,end],.0106,'fitting',3)
    add(f'BPA {bi}/{label} cable leader',[end,A if label=='proximal' else B],.001,'cable',3)
  add('Central Bowden housing',housing,.0025,'bowden',48)
 else:
  # The pressure vessel is straight. Only its output tendon follows the
  # remaining survey guides. Never fold the sleeve around the foot.
  end=L+2*fit;C=P[0]+u*end
  calf=any('calcn' in f for f in r['frames'])
  if calf:
   guide=next((i for i,b in enumerate(r['bodies']) if b=='foot_r'),len(P)-1)
   C=P[guide];u=unit(C-P[0]);O=C-u*end
  else:O=P[0]
  for bi,offset in enumerate(r['offsets'],1):
   off=np.asarray(W[r['offset_body']])[:3,:3]@offset
   off=off-u*np.dot(off,u)
   add(f'BPA {bi}/sleeve',[O+u*fit+off,O+u*(fit+L)+off],D/2,'sleeve',48)
   add(f'BPA {bi}/origin fitting',[O+off,O+u*fit+off],.0106,'fitting',3)
   add(f'BPA {bi}/end fitting',[O+u*(fit+L)+off,C+off],.0106,'fitting',3)
   if calf:add(f'BPA {bi}/origin registration lead',[P[0],O+off],.001,'cable',3)
   add(f'BPA {bi}/collector lead',[C+off,C],.001,'cable',3)
  cable=np.vstack([C,P[guide+1:] if calf else P[1:]])
  if len(cable)<2:cable=np.vstack([C,P[-1]+[0,1e-6,0]])
  add('Shared output tendon',cable,.001,'cable',48)
 # Pinned origin/insertion markers retain the actual recorded/registered
 # attachment coordinate; bundle offsets never move the shared terminal pin.
 for label,O in [('Origin pin',P[0]),('Insertion pin',P[-1])]:
  add(label,[O-np.array([0,0,.004]),O+np.array([0,0,.004])],.004,'pin',3)
 available=float(np.linalg.norm(forward))
 conflict=bool(not r['split_bowden'] and L+2*fit>available+1e-6)
 return out,dict(length_m=L,diameter_m=D,calibration_valid=bool(.75*r['rest']-1e-9<=L<=r['rest']+1e-9),insertion_m=P[-1].tolist(),split_bowden=r['split_bowden'],first_straight_span_m=available,straight_span_conflict=conflict)
