function primitives = fr_sleevePrimitives(C, cs, design, kin, mech, meta)
% Explicit, full-length straight vessels. Placement is a pinned route-tangent
% candidate, not permission to replace a vessel with a zero-length polyline.
N=mech.nG;primitives=cell(1,N);nb=cs.seed.nBPA;
restD=fr_diameter(C,cs.diaStr,0)/1000;
clock=0;
if isfield(cs,'bundleClock_deg'); clock=cs.bundleClock_deg; end
rotate=[cosd(clock) sind(clock);-sind(clock) cosd(clock)];
if isfield(cs,'bundlePacking') && strcmp(cs.bundlePacking,'legacyRing')
    channels=ring(nb)*rotate;
else
    channels=fr_parallelRows(nb,2)*rotate;
end
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
 if isfield(cs,'bpaBendArc_m')
  [parts]=curvedSleeves(P(1:cs.distalSleeveGuideRow,:),L,rad,fit,channels,restD,cs.bundleSpread);
  primitives{i}=parts;
  continue
 end
 split=strcmp(cs.spanKind,'knee')&&design.seriesSplit<.999;
 if split
  origins=[P(1,:);P(end,:)];dirs=[u;v];lengths=L*[design.seriesSplit,1-design.seriesSplit];
 else
  origins=P(1,:);dirs=u;lengths=L;
 end
 parts=struct('p1',{},'p2',{},'r',{},'branch',{},'section',{},'kind',{});
 for section=1:numel(lengths)
  spread=1;if isfield(cs,'bundleSpread');spread=cs.bundleSpread;end
  d=dirs(section,:);perp=null(d);offs=channels*restD/2*spread;
  for b=1:nb
   offset=(perp*offs(b,:)')';a=origins(section,:)+offset+d*fit;
   parts(end+1)=struct('p1',a,'p2',a+d*lengths(section),'r',rad,'branch',b,'section',section,'kind','bpa_sleeve');
  end
 end
 primitives{i}=parts;
end
end

function parts=curvedSleeves(P,L,r,fit,channels,restD,spread)
s=[0;cumsum(vecnorm(diff(P),2,2))];
[s,ids]=unique(s,'stable'); P=P(ids,:);
t=unique([fit;s(s>fit & s<fit+L);fit+L]);
track=interp1(s,P,t,'linear','extrap');
axis=(P(end,:)-P(1,:)); axis=axis/norm(axis);
plane=null(axis); offsets=channels*restD/2*spread;
parts=struct('p1',{},'p2',{},'r',{},'branch',{},'section',{},'kind',{});
for b=1:size(channels,1)
    offset=(plane*offsets(b,:)')';
    for j=1:size(track,1)-1
        parts(end+1)=struct('p1',track(j,:)+offset,'p2',track(j+1,:)+offset, ...
            'r',r,'branch',b,'section',1,'kind','bpa_sleeve'); %#ok<AGROW>
    end
end
end

function off=ring(n)
if n==1;off=[0,0];elseif n==2;off=[0,1;0,-1];else;t=(0:n-1)'*2*pi/n;off=[cos(t),sin(t)]/sin(pi/n);end
end
