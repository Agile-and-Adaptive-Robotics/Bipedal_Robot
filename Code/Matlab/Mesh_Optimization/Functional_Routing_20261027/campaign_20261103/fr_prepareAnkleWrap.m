function chosen=fr_prepareAnkleWrap(C,chosen,kin)
% A guide between the two BPA mounting points deflects the vessel outward.
cs=chosen.caseDef; d=chosen.design;
g=cs.distalSleeveGuideRow;
p=mean(cs.seedLoc([1,g],:),1);
if nargin>=3
    [actual,~,ids]=fr_locForPass(C,cs,d,kin,1,true);
    p=mean(actual(ismember(ids,[1,g]),:),1);
end
cs=fr_addViaGuide(cs,g-1,p);
cs.bpaWrapControlRow=g;
cs.distalSleeveGuideRow=g+1;
cs.shankStopRow=cs.shankStopRow+1;
cs.movable(end+1)=struct('row',g,'cloud','Tibia','radius_mm',45, ...
    'patch',fr_patch(C,'Tibia','tibia_r',p,45));
d.uv(end+1,:)=[0 0 C.standoff_mm(2)];
cs.lb=[cs.lb(1:end-3);-1;-1;C.standoff_mm(1);cs.lb(end-2:end)];
cs.ub=[cs.ub(1:end-3);1;1;C.standoff_mm(2);cs.ub(end-2:end)];
cs.nX=cs.nX+3;
cs.xSeed=fr_packX(cs,d);
chosen.caseDef=cs; chosen.design=d;
end
