"""Orthogonal synapse routes with explicit neuron-interior collision checks."""
import heapq,itertools
import numpy as np
from matplotlib.patches import FancyArrowPatch,Circle
from matplotlib.path import Path

def intersects(a,b,center,size):
 x,y=center;w,h=size;left,right,bottom,top=x-w,x+w,y-h,y+h;eps=1e-7
 if abs(a[0]-b[0])<eps:return left+eps<a[0]<right-eps and min(max(a[1],b[1]),top)-max(min(a[1],b[1]),bottom)>eps
 if abs(a[1]-b[1])<eps:return bottom+eps<a[1]<top-eps and min(max(a[0],b[0]),right)-max(min(a[0],b[0]),left)>eps
 raise ValueError('Route segment is not orthogonal')

def render_routes(ax,bounds,pending):
 step=.06;pad=.12;xs=np.arange(min(c[0]-s[0]for c,s in bounds.items())-.55,max(c[0]+s[0]for c,s in bounds.items())+.6,step);ys=np.arange(min(c[1]-s[1]for c,s in bounds.items())-.55,max(c[1]+s[1]for c,s in bounds.items())+.6,step)
 blocked=np.zeros((len(xs),len(ys)),dtype=bool);used=np.zeros_like(blocked,dtype=float)
 for (x,y),(w,h)in bounds.items():blocked|=(abs(xs[:,None]-x)<w+pad)&(abs(ys[None,:]-y)<h+pad)
 def port(point,other):
  point=np.asarray(point,dtype=float);other=np.asarray(other,dtype=float);owners=[]
  for center,size in bounds.items():
   delta=abs(point-center);w,h=size
   if delta[0]<=w+.13 and delta[1]<=h+.13:owners.append((float(np.linalg.norm(delta)),center,size))
  if not owners:return point,None
  _,center,(w,h)=min(owners);x,y=center;direction=other-center
  if abs(direction[0])/w>abs(direction[1])/h:
   sign=1 if direction[0]>=0 else-1;return np.array([x+sign*w,np.clip(point[1],y-h*.8,y+h*.8)]),(sign,0)
  sign=1 if direction[1]>=0 else-1;return np.array([np.clip(point[0],x-w*.8,x+w*.8),y+sign*h]),(0,sign)
 def anchor(p,normal):
  p=p.copy()
  if normal is not None:p+=np.asarray(normal)*(pad+step*2)
  i=int(np.argmin(abs(xs-p[0])));j=int(np.argmin(abs(ys-p[1])))
  if blocked[i,j]:raise ValueError('No clear synapse port')
  return i,j
 records=[];counter=itertools.count()
 for source,target,exc,color,rad,lw in pending:
  start,ns=port(source,target);end,ne=port(target,source);s=anchor(start,ns);goal=anchor(end,ne)
  initial=(*s,-1);queue=[(0.,next(counter),initial)];dist={initial:0.};parent={};terminal=None
  while queue:
   _,_,state=heapq.heappop(queue);i,j,direction=state
   if(i,j)==goal:terminal=state;break
   for d,(di,dj)in enumerate(((1,0),(-1,0),(0,1),(0,-1))):
    ni,nj=i+di,j+dj
    if not(0<=ni<len(xs)and 0<=nj<len(ys))or blocked[ni,nj]:continue
    key=(ni,nj,d);cost=dist[state]+1.+(3. if direction not in(-1,d)else 0.)+.5*used[ni,nj]
    if cost>=dist.get(key,float('inf')):continue
    dist[key]=cost;parent[key]=state;heuristic=abs(ni-goal[0])+abs(nj-goal[1]);heapq.heappush(queue,(cost+heuristic,next(counter),key))
  if terminal is None:raise ValueError('No collision-free synapse route')
  chain=[];node=terminal
  while True:
   i,j,_=node;chain.append((xs[i],ys[j]));used[i,j]+=1.
   if node==initial:break
   node=parent[node]
  chain.reverse();first=np.array(chain[0]);last=np.array(chain[-1]);points=[start]
  points.append(np.array([start[0],first[1]])if ns is not None and ns[1]else np.array([first[0],start[1]]));points+=list(map(np.array,chain));points.append(np.array([end[0],last[1]])if ne is not None and ne[1]else np.array([last[0],end[1]]));points.append(end)
  compact=[]
  for point in points:
   if compact and np.linalg.norm(point-compact[-1])<1e-9:continue
   while len(compact)>1 and abs(np.cross(compact[-1]-compact[-2],point-compact[-1]))<1e-9:compact.pop()
   compact.append(point)
  for a,b in zip(compact[:-1],compact[1:]):
   for center,size in bounds.items():
    if intersects(a,b,center,size):raise ValueError('Synapse crosses a neuron interior')
  path=Path(np.array(compact),[Path.MOVETO]+[Path.LINETO]*(len(compact)-1));ax.add_patch(FancyArrowPatch(path=path,arrowstyle='-|>'if exc else'-',mutation_scale=11,color=color,lw=lw,zorder=2))
  if not exc:ax.add_patch(Circle(end,.055,color=color,zorder=6))
  records.append(dict(source=list(source),target=list(target),vertices=[p.tolist()for p in compact],neuron_intersections=0))
 ax.set_xlim(min(ax.get_xlim()[0],xs[0]),max(ax.get_xlim()[1],xs[-1]));ax.set_ylim(min(ax.get_ylim()[0],ys[0]),max(ax.get_ylim()[1],ys[-1]))
 return records
