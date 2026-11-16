function fr_tryEDLGuideArrow(outputDir)
% Move the lower tibial guide in the directions marked on the EDL overlay.
C=fr_setup();
C.exactSegmentDistance=true;
C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
branch='compliant';
name='ext_dig_compliant';
stage=fullfile(C.camp,'run_20261105_distal_stages');
prior=fullfile(stage,'ankle_shaft_seed_600',name);
S=load([prior '.mat'],'result');
old=S.result.chosen;
S=load(fullfile(stage,'ankle_mounts_600',[name '.mat']),'result');
base=S.result.chosen;
S=load(fullfile(stage,'ankle_mounts_600','ext_hal_compliant.mat'),'result');
inner=S.result.chosen;
if ~isfolder(outputDir); mkdir(outputDir); end
checkpoint=fullfile(outputDir,[name '.mat']);
reportFile=fullfile(outputDir,[name '.json']);
if isfile(checkpoint) && isfile(reportFile)
    S=load(checkpoint,'result');
    verifyDense(C,kin,S.result.chosen,old,inner,branch,outputDir);
    return
end

base=fr_prepareAnkleWrap(C,base,kin);
base.design.uv(1:end-1,:)=old.design.uv;
base.design.L0=old.design.L0;
base.design.tendon=old.design.tendon;
cs=base.caseDef;
g=find([cs.movable.row]==cs.distalSleeveGuideRow,1);
w=find([cs.movable.row]==cs.bpaWrapControlRow,1);
[p,~,rows]=fr_locForPass(C,old.caseDef,old.design,kin,1,true);
A=p(rows==1,:);
first=p(rows==cs.bpaWrapControlRow,:);
last=p(rows==old.caseDef.bpaBendEndRow,:);
B=p(rows==old.caseDef.distalSleeveGuideRow,:);
u=(first-A)/norm(first-A);
v=(B-last)/norm(B-last);
t=[u(:),v(:)]\(B-A)';
corner=A+t(1)*u;
assert(norm(corner-(B-t(2)*v))<1e-6,'Could not recover previous wrap corner');
base.design.uv(w,:)=projectTarget(cs.movable(w).patch,corner, ...
    base.design.uv(w,:),C.standoff_mm,[1 1 1]);
wp=cs.movable(w).patch;
wrapUV=num2cell(base.design.uv(w,:));
assert(norm(wp.project(wrapUV{:})-corner)<1e-3,'Saved wrap corner is outside the rebuilt patch');

% Specify displacements in pelvis axes at home, then solve in tibia axes.
T=fr_chainT(kin,kin.qDefault,'tibia_r','pelvis');
R=T(1:3,1:3);
targets=[5 5; 10 10; 15 15; 20 20; 30 30; 10 20; 20 10; 10 30; 30 10];
gp=cs.movable(g).patch;
[ehl,innerDetail]=fr_evalCoupledAnkle(C,kin,inner.caseDef,inner.design,branch,11);
target=readtable([prior '_coupled_torque.csv']);
[oldE,oldD]=fr_evalCoupledAnkle(C,kin,old.caseDef,old.design,branch,11);
assert(max(abs(target.ankle_deg-oldE.ankle_deg))<1e-9 && ...
    max(abs(target.mtp_angle_r_deg-oldE.second_deg))<1e-9,'Pose order changed');
oldLayer=layer(C,kin,old,oldE,oldD,inner,ehl,innerDetail);
oldFoot=footRows(oldD);
trials=struct('chosen',{},'coupled',{},'layer',{},'shift_mm',{},'stats',{});
for j=1:size(targets,1)
    shift=[targets(j,1),0,targets(j,2)]/1000;
    d=base.design;
    d.uv(g,:)=projectTarget(gp,B+shift*R,d.uv(g,:),C.standoff_mm,[1 3 1],R);
    uv=num2cell(d.uv(g,:));
    moved=(gp.project(uv{:})-B)*R'*1000;
    if min(moved([1 3]))<1 || abs(moved(2))>5; continue; end
    if ~isempty(trials) && min(vecnorm(vertcat(trials.shift_mm)-moved,2,2))<.25
        continue
    end
    [candidate,d]=fr_materializeAnkleWrap(C,kin,cs,d);
    G=fr_ankleCoupledSetup(C,kin,candidate,11);
    d=fr_refitMountLengths(C,kin,candidate,d,branch,G);
    [e,detail]=fr_evalCoupledAnkle(C,kin,candidate,d,branch,11,G);
    assert(isequal(oldFoot,footRows(detail)),'Foot coordinates moved');
    assert(isequal(oldD.loc(old.caseDef.shankStopRow,:), ...
        detail.loc(candidate.shankStopRow,:)),'Shank stop moved');
    c=struct('caseDef',candidate,'design',d,'reeving',old.reeving);
    L=layer(C,kin,c,e,detail,inner,ehl,innerDetail);
    stats=metrics(e,L,target);
    e.localPass=e.localPass && L.layerPass;
    trials(end+1)=struct('chosen',c,'coupled',e,'layer',L, ... %#ok<AGROW>
        'shift_mm',moved,'stats',stats);
    fprintf('Guide shift [%+.1f %+.1f %+.1f] mm: bone %+.1f mm, fittings %+.1f mm, guide %.1f deg, no force %.1f%%\n', ...
        moved,1000*e.boneMargin_m,1000*L.rigidMargin_m,e.guideTurn_deg,100*mean(e.noForce));
end
assert(~isempty(trials),'The patch produced no anterior/lateral guide moves');
oldStats=metrics(oldE,oldLayer,target);
S=vertcat(trials.stats);
% Prefer feasible routes, then less constraint violation and a gentler guide.
pass=arrayfun(@(trial) trial.coupled.localPass,trials)';
[~,idx]=sortrows([~pass,S(:,1),S(:,2)]);
best=trials(idx(1));
record=jsondecode(fileread([prior '.json']));
record=rmfield(record,{'surrogateEvaluations','initialScore','finalScore'});
record.coupled=best.coupled;
record.layer=best.layer;
record.rest_mm=1000*best.chosen.design.L0;
record.tendon_mm=1000*best.chosen.design.tendon;
record.foundLengthFit=~any(best.coupled.noForce);
record.mountRows=[1 best.chosen.caseDef.distalSleeveGuideRow best.chosen.caseDef.shankStopRow];
record.guideShift_mm=best.shift_mm;
record.baselineFolder='ankle_shaft_seed_600';
record.previous=struct('guideTurn_deg',oldE.guideTurn_deg,'stopTurn_deg',oldE.stopTurn_deg, ...
    'boneMargin_m',oldE.boneMargin_m,'orderingMargin_m',oldLayer.orderingMargin_m, ...
    'rigidMargin_m',oldLayer.rigidMargin_m,'stats',oldStats);
dense=fr_densePilotCase(C,kin,best.chosen.caseDef);
[~,rep]=fr_objective(C,dense,kin,fr_packX(dense,best.chosen.design),branch, ...
    struct('quiet',true,'envelopes',{{}}));
record.localConstraintsPass=rep.feasible;
record.fractionMet=rep.fractionPosesMet;
record.constraints=rep.cinfo;
record.robot_Nm=rep.tauPrimary;
record.human_Nm=rep.tauTarget;
result=struct('chosen',best.chosen,'dense',rep);
save(checkpoint,'result','record','trials','oldStats','-v7');
fid=fopen(reportFile,'w');
assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(record));
fclose(fid);
fr_refreshCoupledRoutePlots(outputDir);
verifyDense(C,kin,best.chosen,old,inner,branch,outputDir);
end

function verifyDense(C,kin,candidate,old,inner,branch,outputDir)
file=fullfile(outputDir,'dense_guide_verification.mat');
if isfile(file); return; end
n=21;
[ehl,innerDetail]=fr_evalCoupledAnkle(C,kin,inner.caseDef,inner.design,branch,n);
cases={old,candidate};
checks=cell(1,numel(cases));
for i=1:numel(cases)
    c=cases{i};
    [e,d]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,branch,n);
    L=layer(C,kin,c,e,d,inner,ehl,innerDetail);
    witness=d.bone;
    witness.ankle_deg=e.ankle_deg(witness.poseIndex);
    witness.mtp_deg=e.second_deg(witness.poseIndex);
    checks{i}=struct('coupled',e,'layer',L,'boneWitness',witness);
    fprintf('Dense case %d: bone %+.3f mm at %.1f / %.1f deg, fitting %+.3f mm, no force %.1f%%, pass %d\n', ...
        i,1000*e.boneMargin_m,witness.ankle_deg,witness.mtp_deg, ...
        1000*L.rigidMargin_m,100*mean(e.noForce),e.localPass && L.layerPass);
end
save(file,'checks','-v7');
end

function uv=projectTarget(patch,point,start,bounds,weights,R)
if nargin<6; R=eye(3); end
lb=[-1 -1 bounds(1)]; ub=[1 1 bounds(2)];
clamp=@(x) max(lb,min(ub,x));
place=@(x) patch.project(x(1),x(2),x(3));
cost=@(x) sum((((place(clamp(x))-point)*R').*weights).^2) + ...
    1e-4*sum(((x-clamp(x))./[1 1 30]).^2);
[u,v,s]=ndgrid(-1:.25:1,-1:.25:1,linspace(bounds(1),bounds(2),5));
seeds=[start;u(:),v(:),s(:)];
value=cellfun(cost,num2cell(seeds,2));
[~,order]=sort(value);
uv=start;
for i=order(1:4)'
    x=fminsearch(cost,seeds(i,:),optimset('Display','off','TolX',1e-6));
    x=clamp(x);
    if cost(x)<cost(uv); uv=x; end
end
end

function p=footRows(detail)
p=detail.loc(ismember(detail.owning,{'talus_r','calcn_r','toes_r'}),:);
end

function L=layer(C,kin,c,e,detail,inner,ehl,innerDetail)
assert(isequal(e.ankle_deg,ehl.ankle_deg) && isequal(e.second_deg,ehl.second_deg));
L=fr_compareAnkleLayers(C,kin,e.poses,detail.sleeves,innerDetail.sleeves, ...
    c.caseDef.fittingLength,inner.caseDef.fittingLength);
L.ehlLocalPass=ehl.localPass;
L.layerPass=L.overlapPairs>0 && min([L.orderingMargin_m,L.rigidMargin_m])>=-1e-4;
end

function s=metrics(e,L,target)
ma=e.ankle_Nm.*sign(target.human_ankle_Nm)-abs(target.human_ankle_Nm);
mb=e.second_Nm.*sign(target.human_second_Nm)-abs(target.human_second_Nm);
worst=min([e.boneMargin_m,e.leadMargin_m,e.shaftCableMargin_m,e.footMargin_m, ...
    e.distalGuideMargin_m,e.wrapMargin_m,e.bundleMargin_m,L.orderingMargin_m,L.rigidMargin_m]);
s=[1000*max(0,-worst)+1000*mean(e.noForce),e.guideTurn_deg, ...
    mean(ma>=-1e-6 & ~e.noForce),mean(mb>=-1e-6 & ~e.noForce),mean(e.noForce)];
end
