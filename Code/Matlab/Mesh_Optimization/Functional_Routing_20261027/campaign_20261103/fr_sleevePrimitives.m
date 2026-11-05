function primitives = fr_sleevePrimitives(C, cs, design, kin, mech, meta)
% Explicit, full-length straight vessels. Placement is a pinned route-tangent
% candidate, not permission to replace a vessel with a zero-length polyline.
N=mech.nG;primitives=cell(1,N);nb=cs.seed.nBPA;
restD=fr_diameter(C,cs.diaStr,0)/1000;
for i=1:N
 q=meta.qg.qList{i};P=zeros(size(meta.loc));
 for j=1:size(P,1)
  T=fr_bodyT(kin,q,meta.rowBodies{j});P(j,:)=(T(1:3,1:3)*meta.loc(j,:)'+T(1:3,4))';
 end
 first=find(vecnorm(P-P(1,:),2,2)>1e-9,1);
 last=find(vecnorm(P-P(end,:),2,2)>1e-9,1,'last');
 assert(~isempty(first)&&~isempty(last),'fr_sleevePrimitives:CoincidentRoute','No sleeve axis for %s',cs.name);
 u=(P(first,:)-P(1,:));u=u/norm(u);
 v=P(last,:)-P(end,:);v=v/norm(v);
 L=mech.sleeveLength_m(i);rad=fr_diameter(C,cs.diaStr,mech.strain_p(i))/2000;
 fit=C.fittingLength;if isfield(cs,'fittingLength');fit=cs.fittingLength;end
 split=strcmp(cs.spanKind,'knee')&&design.seriesSplit<.999;
 if split
  origins=[P(1,:);P(end,:)];dirs=[u;v];lengths=L*[design.seriesSplit,1-design.seriesSplit];
 else
  origins=P(1,:);dirs=u;lengths=L;
 end
 parts=struct('p1',{},'p2',{},'r',{},'branch',{},'section',{},'kind',{});
 for section=1:numel(lengths)
  spread=1;if isfield(cs,'bundleSpread');spread=cs.bundleSpread;end
  d=dirs(section,:);perp=null(d);offs=ring(nb)*restD/2*spread;
  for b=1:nb
   offset=(perp*offs(b,:)')';a=origins(section,:)+offset+d*fit;
   parts(end+1)=struct('p1',a,'p2',a+d*lengths(section),'r',rad,'branch',b,'section',section,'kind','bpa_sleeve');
  end
 end
 primitives{i}=parts;
end
end
function off=ring(n)
if n==1;off=[0,0];elseif n==2;off=[0,1;0,-1];else;t=(0:n-1)'*2*pi/n;off=[cos(t),sin(t)]/sin(pi/n);end
end
