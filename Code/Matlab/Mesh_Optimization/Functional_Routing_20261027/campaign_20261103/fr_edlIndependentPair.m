function fr_edlIndependentPair(folder,priorDense)
% Two identical direct-pull EDL BPAs with separate supports and cable paths.
C=fr_setup();
C.exactSegmentDistance=true;
C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
old=fullfile(C.camp,'run_20261105_distal_stages');
name='ext_dig_compliant';
seedFolder=fullfile(old,'ankle_guide_arrow');
innerFolder=fullfile(old,'ankle_distal_sizing_pilot');
S=load(fullfile(seedFolder,[name '.mat']),'result');
c=S.result.chosen;
S=load(fullfile(innerFolder,'ext_hal_compliant.mat'),'result');
inner=S.result.chosen;
assert(c.caseDef.seed.nBPA==1 && c.caseDef.diaNum==20 && inner.caseDef.diaNum==20);
assert(max(c.caseDef.txConfig.gain)==1);
if ~isfolder(folder); mkdir(folder); end
n=4;
[base,baseD]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,'compliant',n);
common=c.caseDef; common.seed.nBPA=2;
[same,~]=fr_evalCoupledAnkle(C,kin,common,c.design,'compliant',n);
assert(isequal(same.ankle_Nm,2*base.ankle_Nm) && ...
    isequal(same.second_Nm,2*base.second_Nm) && ...
    isequal(same.tendonForce_N,base.tendonForce_N),'Independent-path scaling failed');
fprintf('EDL separate mechanics: exact 2x torque, same tension per tendon at shared centerline.\n');
[~,ehl]=fr_evalCoupledAnkle(C,kin,inner.caseDef,inner.design,'compliant',n);
G=fr_ankleCoupledSetup(C,kin,c.caseDef,n);
[points,frames]=fr_locForPass(C,G.cs,c.design,kin,1);
T=fr_chainT(kin,kin.qDefault,'tibia_r','pelvis');
R=T(1:3,1:3);

if nargin<2
file=fullfile(folder,'independent_pair_screen.mat');
if isfile(file)
    S=load(file,'trials','positions'); trials=S.trials; positions=S.positions;
else
    [X,Z]=ndgrid([0,12,24,38,52,66],[-44,-28,-14,0,14,28,44]);
    positions=[X(:),Z(:)];
    trials=cell(size(positions,1),2);
    for side=1:2
        for j=1:size(positions,1)
            p=points;
            shift=[positions(j,1),0,positions(j,2)]/1000;
            p(1:c.caseDef.distalSleeveGuideRow,:)=p(1:c.caseDef.distalSleeveGuideRow,:)+shift*R;
            stop=c.caseDef.shankStopRow;
            stopMove=.36*shift+[0,0,(2*side-3)*.008];
            p(stop,:)=p(stop,:)+stopMove*R;
            for row=find(strcmp(c.caseDef.rowOwning,'toes_r'))
                U=fr_chainT(kin,kin.qDefault,frames{row},'pelvis');
                p(row,:)=p(row,:)+[0,0,(2*side-3)*.004]*U(1:3,1:3);
            end
            t=evaluate(C,kin,c,p,frames,n);
            t.layer=compare(C,kin,t.d,ehl,t.e.poses,c,inner,true);
            t.stopClearance_m=stopClearance(C,t.d,ehl,c,inner);
            t.score=singleScore(t)+200*sum(shift.^2);
            trials{j,side}=t;
        end
        fprintf('EDL branch %d: evaluated %d independent mount shifts.\n',side,size(positions,1));
    end
    save(file,'trials','positions','-v7');
end

% Screen candidate pairs on the same poses, then re-solve dense finalists.
index=cell(1,2);
for side=1:2
    [~,I]=sort(cellfun(@(t)t.score,trials(:,side)));
    index{side}=I(1:min(18,numel(I)));
end
pairs=zeros(0,3);
for i=index{1}'
    for j=index{2}'
        a=trials{i,1}; b=trials{j,2};
        dist=norm(a.d.loc(1,:)-b.d.loc(1,:));
        if dist<.025; continue; end
        L=compare(C,kin,a.d,b.d,a.e.poses,c,c,false);
        S=stopClearance(C,a.d,b.d,c,c);
        score=a.score+b.score+1000*max([0,-L.rigidMargin_m,L.sleeveOverlap_m,-S]);
        pairs(end+1,:)=[score,i,j]; %#ok<AGROW>
    end
end
assert(~isempty(pairs),'No two distinct EDL mount candidates');
pairs=sortrows(pairs,1);
n=21;
[baseline,~]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,'compliant',n);
[ehlE,ehlD]=fr_evalCoupledAnkle(C,kin,inner.caseDef,inner.design,'compliant',n);
best=inf; chosen=[]; denseTrials={};
for k=1:min(3,size(pairs,1))
    a=evaluate(C,kin,c,trials{pairs(k,2),1}.d.loc,frames,n);
    b=evaluate(C,kin,c,trials{pairs(k,3),2}.d.loc,frames,n);
    a.layer=compare(C,kin,a.d,ehlD,a.e.poses,c,inner,true);
    b.layer=compare(C,kin,b.d,ehlD,b.e.poses,c,inner,true);
    a.stopClearance_m=stopClearance(C,a.d,ehlD,c,inner);
    b.stopClearance_m=stopClearance(C,b.d,ehlD,c,inner);
    L=compare(C,kin,a.d,b.d,a.e.poses,c,c,false);
    S=stopClearance(C,a.d,b.d,c,c);
    score=singleScore(a)+singleScore(b)+1000*max([0,-L.rigidMargin_m,L.sleeveOverlap_m,-S]);
    t=struct('a',a,'b',b,'pair',L,'stopClearance_m',S,'score',score);
    denseTrials{end+1}=t; %#ok<AGROW>
    if score<best; best=score; chosen=t; end
end
else
    S=load(priorDense,'denseTrials');
    [~,ehlCoarse]=fr_evalCoupledAnkle(C,kin,inner.caseDef,inner.design,'compliant',5);
    [a,b,stopTrials]=refineStops(C,kin,c,inner,S.denseTrials,frames,R,ehlCoarse);
    save(fullfile(folder,'stop_refinement.mat'),'stopTrials','-v7');
    n=21;
    [baseline,~]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,'compliant',n);
    [ehlE,ehlD]=fr_evalCoupledAnkle(C,kin,inner.caseDef,inner.design,'compliant',n);
    a=evaluate(C,kin,c,a.d.loc,frames,n);
    b=evaluate(C,kin,c,b.d.loc,frames,n);
    a.layer=compare(C,kin,a.d,ehlD,a.e.poses,c,inner,true);
    b.layer=compare(C,kin,b.d,ehlD,b.e.poses,c,inner,true);
    a.stopClearance_m=stopClearance(C,a.d,ehlD,c,inner);
    b.stopClearance_m=stopClearance(C,b.d,ehlD,c,inner);
    L=compare(C,kin,a.d,b.d,a.e.poses,c,c,false);
    S=stopClearance(C,a.d,b.d,c,c);
    wire=wireMargin(kin,c,a,b);
    chosen=struct('a',a,'b',b,'pair',L,'stopClearance_m',S,'wireMargin_m',wire);
    denseTrials={chosen};
end
a=chosen.a; b=chosen.b;
e=a.e;
e.ankle_Nm=a.e.ankle_Nm+b.e.ankle_Nm;
e.second_Nm=a.e.second_Nm+b.e.second_Nm;
e.tendonForce_N=a.e.tendonForce_N+b.e.tendonForce_N;
e.totalTendonForce_N=e.tendonForce_N;
e.slack=a.e.slack & b.e.slack;
e.stroke=a.e.stroke & b.e.stroke;
e.invalid=a.e.invalid & b.e.invalid;
e.noForce=a.e.noForce & b.e.noForce;
e.localPass=a.e.localPass && b.e.localPass && ...
    min([a.layer.rigidMargin_m,b.layer.rigidMargin_m,a.layer.orderingMargin_m, ...
    b.layer.orderingMargin_m,chosen.pair.rigidMargin_m,a.stopClearance_m, ...
    b.stopClearance_m,chosen.stopClearance_m])>=-1e-4 && chosen.pair.sleeveOverlap_m<=1e-9;
if isfield(chosen,'wireMargin_m'); e.localPass=e.localPass && chosen.wireMargin_m>=-1e-4; end
e.boneMargin_m=min(a.e.boneMargin_m,b.e.boneMargin_m);
e.leadMargin_m=min(a.e.leadMargin_m,b.e.leadMargin_m);
e.bundleMargin_m=-chosen.pair.sleeveOverlap_m;
record=jsondecode(fileread(fullfile(seedFolder,[name '.json'])));
record.legacySweep=false;
record.count=2;
record.independentLoadPaths=true;
record.coupled=e;
record.baselineCoupled=baseline;
record.baselineFolder='';
record.foundLengthFit=~any(a.e.noForce | b.e.noForce) && ...
    c.design.tendon>=max(a.e.tendonRequired_m,b.e.tendonRequired_m)-1e-6;
record.localConstraintsPass=e.localPass;
record.ehlLocalPass=ehlE.localPass;
record.assembledLocalPass=e.localPass && ehlE.localPass;
record.layer=struct('rigidMargin_m',min(a.layer.rigidMargin_m,b.layer.rigidMargin_m), ...
    'orderingMargin_m',min(a.layer.orderingMargin_m,b.layer.orderingMargin_m), ...
    'sleeveOverlap_m',max(a.layer.sleeveOverlap_m,b.layer.sleeveOverlap_m));
record.pairClearance=chosen.pair;
record.stopClearance_m=chosen.stopClearance_m;
if isfield(chosen,'wireMargin_m'); record.wireMargin_m=chosen.wireMargin_m; end
record.routeParts={nativeRoute(kin,a.d),nativeRoute(kin,b.d)};
record.points_m=record.routeParts{1}.points_m;
record.bodies=record.routeParts{1}.bodies;
record.innerRoute=nativeRoute(kin,ehlD);
record.branchResults={a.e,b.e};
record.branchBoneWitnesses={a.d.bone,b.d.bone};
record.mountMovement_mm={1000*(a.d.loc([1,c.caseDef.distalSleeveGuideRow,c.caseDef.shankStopRow],:)- ...
    baseD.loc([1,c.caseDef.distalSleeveGuideRow,c.caseDef.shankStopRow],:))*R', ...
    1000*(b.d.loc([1,c.caseDef.distalSleeveGuideRow,c.caseDef.shankStopRow],:)- ...
    baseD.loc([1,c.caseDef.distalSleeveGuideRow,c.caseDef.shankStopRow],:))*R'};
record=rmfield(record,intersect(fieldnames(record),{'guideShift_mm','previous'}));
save(fullfile(folder,'independent_pair_dense.mat'),'chosen','denseTrials','record','c','inner','-v7');
fid=fopen(fullfile(folder,[name '.json']),'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(record)); fclose(fid);
fprintf('EDL pair dense: bone %+.2f mm, EHL fittings %+.2f mm, mutual fittings %+.2f mm, combined no force %.1f%%\n', ...
    1000*e.boneMargin_m,1000*record.layer.rigidMargin_m,1000*chosen.pair.rigidMargin_m,100*mean(e.noForce));
end

function t=evaluate(C,kin,c,p,frames,n)
G=fr_ankleCoupledSetup(C,kin,c.caseDef,n);
G.qg.resolvedRoute=struct('points',p,'frames',{frames});
[e,d]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,'compliant',n,G);
t=struct('e',e,'d',d);
end

function L=compare(C,kin,a,b,poses,outer,inner,ordered)
L=fr_compareAnkleLayers(C,kin,poses,a.sleeves,b.sleeves, ...
    outer.caseDef.fittingLength,inner.caseDef.fittingLength,ordered);
end

function margin=stopClearance(C,a,b,ca,cb)
% Separate small cable-stop envelopes on the shank.
A=a.loc(ca.caseDef.shankStopRow,:);
B=b.loc(cb.caseDef.shankStopRow,:);
diameter=hypot(C.bowden.stopAllowance_mm,C.bowden.housingOD_mm)/1000;
margin=norm(A-B)-diameter-C.clearanceEnvelope_m;
end

function f=singleScore(t)
e=t.e;
margin=[e.boneMargin_m,e.leadMargin_m,e.footMargin_m,e.distalGuideMargin_m, ...
    e.wrapMargin_m,e.shaftCableMargin_m,t.layer.rigidMargin_m, ...
    t.layer.orderingMargin_m,t.stopClearance_m];
f=100*mean(e.noForce)+1000*sum(max(0,-margin))+50*t.layer.sleeveOverlap_m;
end

function [a,b,trials]=refineStops(C,kin,c,inner,seeds,frames,R,ehl)
n=5;
G=fr_ankleCoupledSetup(C,kin,c.caseDef,n);
stop=c.caseDef.shankStopRow;
ankle=fr_chainT(kin,kin.qDefault,'talus_r','tibia_r');
[X,Y,Z]=ndgrid([-20,-10,0,10,20],[-12,0,10],[-32,-24,-16,-8,0,8]);
offsets=[X(:),Y(:),Z(:)]/1000;
trials=cell(numel(seeds),2);
best=inf; a=[]; b=[];
for k=1:numel(seeds)
    paths={seeds{k}.a.d.loc,seeds{k}.b.d.loc};
    for row=find(strcmp(c.caseDef.rowOwning,'calcn_r'))
        if norm(paths{1}(row,:)-paths{2}(row,:))>=.002; continue; end
        U=fr_chainT(kin,kin.qDefault,frames{row},'pelvis');
        for side=1:2
            paths{side}(row,:)=paths{side}(row,:)+[0,0,(2*side-3)*.004]*U(1:3,1:3);
        end
    end
    for side=1:2
        p=paths{side};
        rows=zeros(0,2);
        for j=1:size(offsets,1)
            point=p; point(stop,:)=point(stop,:)+offsets(j,:)*R;
            if point(stop,2)<ankle(2,4)+.003; continue; end
            span=fr_spans(C,kin,point,frames,G.cs.passes{1}.qStructs,'foot');
            deficit=max(0,span.requiredSpan_m-c.design.tendon);
            rows(end+1,:)=[1000*deficit+norm(offsets(j,:)),j]; %#ok<AGROW>
        end
        rows=sortrows(rows,1);
        pool=cell(1,min(10,size(rows,1)));
        for i=1:numel(pool)
            j=rows(i,2);
            point=p; point(stop,:)=point(stop,:)+offsets(j,:)*R;
            t=evaluate(C,kin,c,point,frames,n);
            t.layer=compare(C,kin,t.d,ehl,t.e.poses,c,inner,true);
            t.stopClearance_m=stopClearance(C,t.d,ehl,c,inner);
            t.score=singleScore(t)+1000*max(0,t.e.tendonRequired_m-c.design.tendon);
            pool{i}=t;
        end
        trials{k,side}=pool;
    end
    for A=trials{k,1}
        for B=trials{k,2}
            x=A{1}; y=B{1};
            L=compare(C,kin,x.d,y.d,x.e.poses,c,c,false);
            S=stopClearance(C,x.d,y.d,c,c);
            wire=wireMargin(kin,c,x,y);
            score=x.score+y.score+1000*max([0,-L.rigidMargin_m,L.sleeveOverlap_m,-S,-wire]);
            if score<best; best=score; a=x; b=y; end
        end
    end
end
assert(~isempty(a),'No stop layout above the ankle');
fprintf('Stop refinement: predicted margins to the fixed tendon lengths %+.2f and %+.2f mm.\n', ...
    1000*(c.design.tendon-a.e.tendonRequired_m),1000*(c.design.tendon-b.e.tendonRequired_m));
end

function margin=wireMargin(kin,c,a,b)
% Wire-to-wire clearance from the independent stops through the toe routes.
assert(isequal(a.e.poses,b.e.poses));
rows=c.caseDef.shankStopRow:size(a.d.loc,1);
margin=inf;
for i=1:numel(a.e.poses)
    A=zeros(numel(rows),3); B=A;
    for j=1:numel(rows)
        r=rows(j);
        T=fr_bodyT(kin,a.e.poses{i},a.d.frames{r});
        A(j,:)=(T(1:3,1:3)*a.d.loc(r,:)'+T(1:3,4))';
        B(j,:)=(T(1:3,1:3)*b.d.loc(r,:)'+T(1:3,4))';
    end
    for j=1:size(A,1)-1
        for k=1:size(B,1)-1
            margin=min(margin,fr_segmentDistance(A(j,:),A(j+1,:),B(k,:),B(k+1,:))-.002);
        end
    end
end
end

function r=nativeRoute(kin,d)
p=d.loc;
for j=1:size(p,1)
    T=fr_chainT(kin,kin.qDefault,d.frames{j},d.owning{j});
    p(j,:)=(T(1:3,1:3)*p(j,:)'+T(1:3,4))';
end
r=struct('points_m',p,'bodies',{d.owning});
end
