function [points,geometry] = fr_seedAnkleShaftWrap(C,kin,cs,d)
% Propose routes from bounded bone patches before the coupled surrogate.
% Several shaft stations are paired with a stop near the continuation
% of the incoming cable. This first screen is geometric, at home.
coupled=fr_ankleCoupledSetup(C,kin,cs,2);
[P,frames,ids]=fr_locForPass(C,coupled.cs,d,kin,1,true);
foot=find(ismember(cs.rowOwning,{'talus_r','calcn_r','toes_r'}),1);
row=cs.shankStopRow;
assert(foot==row+1);
T=fr_chainT(kin,kin.qDefault,frames{ids==foot},'tibia_r');
destination=(T(1:3,1:3)*P(ids==foot,:)'+T(1:3,4))';
guideIndex=find([cs.movable.row]==cs.distalSleeveGuideRow);
stopIndex=find([cs.movable.row]==row);
wrapIndex=find([cs.movable.row]==cs.bpaWrapControlRow);
guidePatch=cs.movable(guideIndex).patch;
stopPatch=cs.movable(stopIndex).patch;
wrapPatch=cs.movable(wrapIndex).patch;
guides=pool(guidePatch,d.uv(guideIndex,:),C.standoff_mm);
stops=pool(stopPatch,d.uv(stopIndex,:),C.standoff_mm);
wraps=pool(wrapPatch,d.uv(wrapIndex,:),C.standoff_mm);
origin=P(ids==1,:);
guides=guides(guides(:,5)<origin(2)-.06,:);
points=zeros(0,cs.nX); geometry=zeros(0,3);
for j=1:size(guides,1)
    g=guides(j,4:6);
    candidates=stops(stops(:,5)<g(2)-.012 & stops(:,5)>destination(2),:);
    if isempty(candidates); continue; end
    v=destination-g;
    t=max(0,min(1,((candidates(:,4:6)-g)*v')/max(dot(v,v),eps)));
    offset=vecnorm(candidates(:,4:6)-g-t.*v,2,2);
    [~,rank]=sort(offset);
    shaft=wraps(wraps(:,5)>g(2)+.02 & wraps(:,5)<origin(2)-.02,:);
    if isempty(shaft); continue; end
    bend=shaft(:,4:6)-origin;
    u=g-origin;
    dist=vecnorm(bend-((bend*u')/dot(u,u)).*u,2,2);
    [~,pick]=sort(dist);
    pick=unique(pick(round(linspace(1,numel(pick),3))));
    for k=rank(1:min(2,numel(rank)))'
        s=candidates(k,4:6);
        x=[g;s;destination];
        margin=inf;
        for q=1:2
            t0=linspace(.12,.88,9)';
            segment=x(q,:)+t0.*(x(q+1,:)-x(q,:));
            margin=min(margin,min(fr_surfaceDistance(C,'tibia_r',segment))-.002);
        end
        for w=pick'
            v0=d;
            v0.uv(guideIndex,:)=guides(j,1:3);
            v0.uv(stopIndex,:)=candidates(k,1:3);
            v0.uv(wrapIndex,:)=shaft(w,1:3);
            p0=shaft(w,4:6);
            angles=[angle(p0-origin,g-p0),angle(g-p0,s-g),angle(s-g,destination-s)];
            % Keep clearance first, then favor the gentler transitions.
            cost=1e5*max(0,-margin)+sum(angles.^2);
            points(end+1,:)=fr_packX(cs,v0)'; %#ok<AGROW>
            geometry(end+1,:)=[cost,margin,max(angles)]; %#ok<AGROW>
        end
    end
end
[~,order]=sort(geometry(:,1));
[~,distinct]=unique(round(points(order,:),5),'rows','stable');
idx=order(distinct(1:min(24,numel(distinct))));
points=points(idx,:);
geometry=geometry(idx,:);
end

function p=pool(patch,seed,standoff)
[U,V,S]=ndgrid([-.8,0,.8],[-.8,0,.8],linspace(standoff(1),standoff(2),3));
use=U.^2+V.^2<=1;
p=[U(use),V(use),S(use);seed];
p(:,4:6)=0;
for i=1:size(p,1)
    p(i,4:6)=patch.project(p(i,1),p(i,2),p(i,3));
end
end

function a=angle(u,v)
a=acos(max(-1,min(1,dot(u,v)/max(norm(u)*norm(v),eps))));
end
