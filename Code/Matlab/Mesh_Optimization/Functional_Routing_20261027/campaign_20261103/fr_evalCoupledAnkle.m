function [E,detail] = fr_evalCoupledAnkle(C,kin,caseDef,design,branch,n,G)
if nargin<7; G=fr_ankleCoupledSetup(C,kin,caseDef,n); end
cs=G.cs;
[m,meta]=fr_mechBuild(C,cs,design,kin,branch,G.qg);
E.ankle_Nm=fr_coordinateTorque(cs,kin,m,1,1,G.indices);
cs.passes{1}.sweptCoord=G.second;
E.second_Nm=fr_coordinateTorque(cs,kin,m,1,2,G.indices);
E.slack=logical(m.cell.PulleySlack(:));
E.stroke=logical(m.cell.PulleyInfeasible(:)) | ...
    m.cell.sContraction(:)/m.rest > 1-m.kmax/m.rest+1e-9;
E.invalid=~isfinite(E.ankle_Nm)|~isfinite(E.second_Nm);
E.tendonForce_N=m.cell.Ftendon(:);
E.totalTendonForce_N=m.forceCopies*E.tendonForce_N;
E.noForce=E.slack|E.stroke|E.invalid|~isfinite(E.tendonForce_N)|E.tendonForce_N<=1e-6;
E.ankle_Nm(E.slack|E.stroke|E.invalid)=0;
E.second_Nm(E.slack|E.stroke|E.invalid)=0;
E.ankle_deg=G.ankle_deg; E.second_deg=G.second_deg; E.second=G.second;
E.poses=G.cs.passes{1}.qStructs;

% Evaluate the same solved sleeves and free lead at each combined pose.
poses=E.poses;
meta.qg.qList=poses;
m.nG=numel(poses);
m.strain_p=m.cell.sContraction(:)/m.rest;
m.sleeveLength_m=m.rest-m.cell.sContraction(:);
sleeves=fr_sleevePrimitives(C,cs,design,kin,m,meta);
obs=fr_obstacles(C,cs,meta.loc,meta.rowBodies,kin,poses,{},sleeves);
E.boneMargin_m=obs.boneMargin;
if nargout>1
    [E.leadMargin_m,lead]=fr_tibialLeadMargin(C,kin,meta.loc,meta.rowBodies, ...
        cs.rowOwning,poses,.002);
    detail=struct('bone',obs.boneWitness,'lead',lead,'sleeves',{sleeves}, ...
        'owning',{cs.rowOwning},'loc',meta.loc,'frames',{meta.rowBodies}, ...
        'strain',m.strain_p,'length_m',m.sleeveLength_m);
else
    E.leadMargin_m=fr_tibialLeadMargin(C,kin,meta.loc,meta.rowBodies,cs.rowOwning,poses,.002);
end
E.tendonRequired_m=fr_spans(C,kin,meta.loc,meta.rowBodies,poses,'foot').requiredSpan_m;
E.footMargin_m=inf; E.bundleMargin_m=inf; E.bundleHardwareMargin_m=inf;
E.distalGuideMargin_m=inf;
if isfield(cs,'distalSleeveGuideRow')
    guide=sum(vecnorm(diff(meta.loc(1:cs.distalSleeveGuideRow,:)),2,2));
    E.distalGuideMargin_m=guide-2*cs.fittingLength-max(m.sleeveLength_m);
end
E.wrapMargin_m=inf;
E.shaftCableMargin_m=inf;
if isfield(cs,'bpaBendArc_m') && cs.bpaBendArc_m>0
    bendEnd=sum(vecnorm(diff(meta.loc(1:cs.bpaBendEndRow,:)),2,2));
    bendStart=bendEnd-cs.bpaBendArc_m;
    E.wrapMargin_m=min([cs.bpaBendRadius_m-cs.bpaBendRadiusTarget_m, ...
        min(m.sleeveLength_m)+cs.fittingLength-bendEnd, ...
        bendStart-cs.fittingLength]);
    E.bpaBendRadius_m=cs.bpaBendRadius_m;
    E.bpaBendArc_m=cs.bpaBendArc_m;
    E.bpaBendLossMax_m=max(m.cell.deltaRoute(:));
end
if isfield(cs,'distalSleeveGuideRow') && isfield(cs,'shankStopRow') && ...
        (isfield(cs,'assemblySizing') || isfield(cs,'bpaBendArc_m'))
    L=fr_ankleShaftLead(C,kin,cs,meta.loc,meta.rowBodies,poses);
    E.shaftCableMargin_m=L.shaftCableMargin_m;
    E.guideTurn_deg=L.guideTurn_deg;
    E.stopTurn_deg=L.stopTurn_deg;
end
for i=1:numel(poses)
    T=fr_bodyT(kin,poses{i},'tibia_r');
    cap=sleeves{i};
    for a=1:numel(cap)
        x=cap(a);
        pts=T(1:3,1:3)'*([x.p1;x.p2]'-T(1:3,4));
        E.footMargin_m=min(E.footMargin_m,min(pts(2,:))+.43-x.r-cs.fittingLength-.005);
        for b=a+1:numel(cap)
            y=cap(b);
            if x.branch==y.branch && x.section==y.section; continue; end
            d=fr_segmentDistance(x.p1,x.p2,y.p1,y.p2)-x.r-y.r;
            E.bundleMargin_m=min(E.bundleMargin_m,d);
        end
    end
    if cs.seed.nBPA>1
        [centers,radii,branches]=fr_fittingSpheres(cap,cs.fittingLength);
        for a=1:size(centers,1)
            for b=find(branches~=branches(a))'
                gap=norm(centers(a,:)-centers(b,:))-radii(a)-radii(b);
                E.bundleHardwareMargin_m=min(E.bundleHardwareMargin_m,gap-C.clearanceEnvelope_m);
            end
            for b=cap([cap.branch]~=branches(a))
                gap=fr_segmentDistance(centers(a,:),centers(a,:),b.p1,b.p2)-radii(a)-b.r;
                E.bundleHardwareMargin_m=min(E.bundleHardwareMargin_m,gap-C.clearanceEnvelope_m);
            end
        end
    end
end
E.localPass=~any(E.noForce) && ...
    min([E.boneMargin_m,E.leadMargin_m,E.footMargin_m,E.distalGuideMargin_m, ...
        E.wrapMargin_m,E.shaftCableMargin_m,E.bundleHardwareMargin_m])>=-1e-4 && ...
    E.bundleMargin_m>=-1e-9 && design.tendon>=E.tendonRequired_m-1e-6;
end
