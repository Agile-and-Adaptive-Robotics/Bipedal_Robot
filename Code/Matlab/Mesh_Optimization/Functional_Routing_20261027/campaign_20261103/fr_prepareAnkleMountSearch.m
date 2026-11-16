function chosen=fr_prepareAnkleMountSearch(C,kin,chosen)
% Three independent tibial points: proximal mount, sleeve-axis guide, stop.
cs=chosen.caseDef; d=chosen.design;
foot=find(ismember(cs.rowOwning,{'talus_r','calcn_r','toes_r'}),1);
assert(foot>=2 && all(strcmp(cs.rowOwning(1:foot-1),'tibia_r')));
if foot==2
    T=fr_chainT(kin,kin.qDefault,cs.rowOwning{foot},'tibia_r');
    p=(T(1:3,1:3)*cs.seedLoc(foot,:)'+T(1:3,4))';
    s=cs.seedLoc(1,:)+.86*(p-cs.seedLoc(1,:));
    cs=fr_addViaGuide(cs,1,s,true);
    foot=foot+1;
end
if foot==3
    % Keep the initial straight lead collinear before projecting the new pin.
    p=.25*cs.seedLoc(1,:)+.75*cs.seedLoc(2,:);
    cs=fr_addViaGuide(cs,1,p);
    foot=foot+1;
end
cs.distalSleeveGuideRow=2;
cs.shankStopRow=foot-1;
for row=[1 cs.distalSleeveGuideRow cs.shankStopRow]
    if any([cs.movable.row]==row); continue; end
    m=struct('row',row,'cloud','Tibia','radius_mm',25, ...
        'patch',fr_patch(C,'Tibia','tibia_r',cs.seedLoc(row,:),25));
    cs.movable(end+1)=m;
    d.uv(end+1,:)=[0 0 mean(C.standoff_mm)];
    cs.lb=[cs.lb(1:end-3);-1;-1;C.standoff_mm(1);cs.lb(end-2:end)];
    cs.ub=[cs.ub(1:end-3);1;1;C.standoff_mm(2);cs.ub(end-2:end)];
    cs.nX=cs.nX+3;
end
cs.xSeed=fr_packX(cs,d);
chosen.caseDef=cs; chosen.design=d;
end
