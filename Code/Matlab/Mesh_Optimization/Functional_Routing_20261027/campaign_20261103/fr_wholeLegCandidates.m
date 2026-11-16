function fr_wholeLegCandidates(out,mode)
% Build independently solved route alternatives on common right-leg poses.
C=fr_setup();
kin=fr_robotKneeContext(C,fr_kinematics(C.osimStock));
C.exactSegmentDistance=true; C.fullFittingEnvelopes=true;
upper=fullfile(C.camp,'run_20261106_remaining20');
ankle=fullfile(C.camp,'run_20261106_ankle_assembled','right_side');
A=load(fullfile(ankle,'assembled_validation.mat'),'routes');
tc=jsondecode(fileread(fullfile(C.targetsDir,'target_contract.json')));
groups={tc.groups.group};
keys={'hip_flexion_r','hip_adduction_r','hip_rotation_r','knee_angle_r', ...
    'ankle_angle_r','subtalar_angle_r','mtp_angle_r','lumbar_extension','lumbar_bending','lumbar_rotation'};
rom=jsondecode(fileread(fullfile(C.repo,'Code','Matlab','Mesh_Optimization','Routing27_v2','actuator_map_27.json'))).rom_table_deg;
limits=zeros(numel(keys),2);
for j=1:numel(keys)
    if isfield(rom,keys{j}); limits(j,:)=rom.(keys{j}); else; limits(j,:)=[-20 20]; end
end
limits(4,:)=[-120,0];
num=24; sweep=5;
if strcmp(mode,'dense'); num=96; sweep=13; end
angles=zeros(1+num+numel(keys)*sweep,numel(keys));
primes=[2 3 5 7 11 13 17 19 23 29];
for j=1:numel(keys)
    angles(1,j)=min(limits(j,2),max(limits(j,1),rad2deg(kin.qDefault.(keys{j}))));
    for i=1:num
        angles(i+1,j)=limits(j,1)+radical(i,primes(j))*diff(limits(j,:));
    end
end
for j=1:numel(keys)
    inds=1+num+(j-1)*sweep+(1:sweep);
    angles(inds,:)=repmat(angles(1,:),sweep,1);
    angles(inds,j)=linspace(limits(j,1),limits(j,2),sweep);
end
poses=cell(1,size(angles,1));
for i=1:numel(poses)
    q=kin.qDefault;
    for j=1:numel(keys); q.(keys{j})=deg2rad(angles(i,j)); end
    poses{i}=q;
end
root=fullfile(out,mode);
if ~isfolder(root); mkdir(root); end
bodies={'pelvis','femur_r','tibia_r','talus_r','calcn_r','toes_r','torso'};
world=struct();
for b=bodies
    transforms=zeros(numel(poses),4,4);
    for i=1:numel(poses); transforms(i,:,:)=fr_bodyT(kin,poses{i},b{1}); end
    world.(b{1})=transforms;
end
write(fullfile(root,'poses.json'),struct('keys',{keys},'angles_deg',angles, ...
    'transforms',world,'reference',C.osimRobot,'anatomy',C.osimStock));
if strcmp(mode,'dense')
    selection=jsondecode(fileread(fullfile(out,'selection.json')));
else
    selection=struct();
end
pool=gcp('nocreate');
if isempty(pool)
    cluster=parcluster('local');
    cluster.JobStorageLocation=fullfile(out,['pack_workers_' mode]);
    if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
    parpool(cluster,min(8,cluster.NumWorkers));
end
branches={'rigid','compliant'};
[j,k]=ndgrid(1:numel(groups),1:2);
jobs=[j(:),k(:)];
parfor i=1:size(jobs,1)
    g=groups{jobs(i,1)}; b=branches{jobs(i,2)};
    file=fullfile(root,[g '_' b '.json']);
    if isfile(file); continue; end
    try
        prepare(C,kin,upper,A.routes,g,b,poses,root,selection);
    catch ME
        fid=fopen(fullfile(root,[g '_' b '_error.txt']),'w');
        fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
        fprintf('FAILED %s %s: %s\n',g,b,ME.message);
    end
end
assert(numel(dir(fullfile(root,'*_rigid.json')))==27 && ...
       numel(dir(fullfile(root,'*_compliant.json')))==27,'Incomplete packaging candidate set');
end

function prepare(C,kin,upper,ankle,g,branch,poses,root,selection)
id=[g '_' branch];
routes=ankle(cellfun(@(x)strcmp(x.group,g),ankle));
isAnkle=~isempty(routes);
targetFile=fullfile(upper,'targets',[g '_search.json']);
target=[];
if isfile(targetFile); target=jsondecode(fileread(targetFile)); end
if ~isAnkle
    file=fullfile(upper,'paired_refit',[id '.mat']);
    if ~isfile(file); file=fullfile(upper,[id '.mat']); end
    S=load(file,'chosen','record');
    r=struct('caseDef',S.chosen.caseDef,'design',S.chosen.design, ...
        'p',S.chosen.caseDef.seedLoc,'channel',1);
    ids=find(r.caseDef.rowUse{1});
    r.p(ids,:)=S.record.points_m;
    routes={r};
end
for j=1:numel(routes)
    routes{j}.caseDef.xi3Mode='perPose';
    if isfield(routes{j}.caseDef.txConfig,'pulleyExitIndex')
        routes{j}.caseDef.txConfig=rmfield(routes{j}.caseDef.txConfig,'pulleyExitIndex');
    end
    routes{j}.dense=fr_densePilotCase(C,kin,routes{j}.caseDef);
end
% Lateral/anterior lane changes and axially staggered fittings.
offsets=[0 0 0;12 0 0;-12 0 0;0 0 12;0 0 -12; ...
    24 0 0;-24 0 0;0 0 24;0 0 -24;12 16 12; ...
    -12 -16 12;12 -16 -12;-12 16 -12;0 20 0;0 -20 0; ...
    20 0 20;-20 0 20;20 0 -20;-20 0 -20];
variants=[offsets,ones(size(offsets,1),1); ...
    0 0 0 .75;0 0 0 1.25;12 0 0 .85;0 0 12 .85];
if isAnkle; variants=variants(1:19,:); end
ids=1:size(variants,1);
if ~isempty(fieldnames(selection)); ids=unique([1,selection.(id)],'stable'); end
options=cell(1,numel(ids)); baseline=[];
for k=1:numel(ids)
    index=ids(k);
    variant=variants(index,:);
    candidate=routes;
    for j=1:numel(candidate)
        r=candidate{j};
        r=relocate(C,kin,r,variant,isAnkle);
        candidate{j}=r;
    end
    result=evaluate(C,kin,candidate,g,branch,poses,isAnkle,target);
    if index==1
        baseline=result;
    elseif any(result.capacity<baseline.capacity-5)
        % Fit only after an attachment change costs measurable torque.
        for j=1:numel(candidate)
            r=candidate{j};
            r.design.tendon=max(.008,r.design.tendon + ...
                routeLength(kin,r.p,r.caseDef.rowOwning)- ...
                routeLength(kin,routes{j}.p,r.caseDef.rowOwning));
            candidate{j}=r;
        end
        alternate=evaluate(C,kin,candidate,g,branch,poses,isAnkle,target);
        if alternate.capacity>result.capacity && alternate.active>=result.active
            result=alternate;
        end
    end
    result.id=index;
    result.offset=variant;
    options{k}=result;
end
write(fullfile(root,[id '.json']),options);
fprintf('%s %s: %d alternatives on %d common poses\n',g,branch,numel(ids),numel(poses));
end

function r=relocate(C,kin,r,variant,isAnkle)
cs=r.caseDef;
if ~isAnkle
    cs.bundleSpread=max(.8,variant(4)*cs.bundleSpread);
end
v=variant(1:3)/1000;
if any(v)
    first=cs.rowOwning{1};
    if isAnkle
        ids=find(strcmp(cs.rowOwning,'tibia_r'));
    else
        ids=intersect([cs.movable.row],find(strcmp(cs.rowOwning,first)),'stable');
        if isempty(ids); ids=1; end
    end
    r.p(ids,:)=r.p(ids,:)+v;
    if isAnkle && any(isfield(cs,{'bpaBendArc_m'})) && isfield(cs,'bpaBendEndRow')
        % A parallel shift preserves the saved curved sleeve's bend.
        assert(all(ismember(1:cs.bpaBendEndRow,ids)));
    end
end
r.caseDef=cs;
end

function D=evaluate(C,kin,routes,g,branch,poses,isAnkle,target)
N=numel(poses); parts=cell(1,N); lines=cell(1,numel(routes));
active=false(N,1); peak=0; minBone=inf; minWire=inf; home=inf;
bay=inf; tendon=inf; footBay=inf;
qHome=kin.qDefault;
for key=fieldnames(qHome)'; qHome.(key{1})=0; end
primary=zeros(31,1); hasForce=false(31,1); targetPrimary=[];
for j=1:numel(routes)
    r=routes{j}; cs=r.caseDef;
    [m,meta,geom]=atPoses(C,kin,r,branch,poses,isAnkle,target);
    vessel=fr_sleevePrimitives(C,geom,r.design,kin,m,meta);
    active=active | transmission(m);
    peak=max(peak,max([0;m.Ftendon(:)],[],'omitnan'));
    dense=r.dense;
    h=dense.target.tau(:);
    [t,a]=primaryCurve(C,kin,r,dense,branch,isAnkle);
    primary=primary+t; hasForce=hasForce|a; targetPrimary=h;
    p=r.p;
    if ~strcmp(cs.spanKind,'none')
        span=fr_spans(C,kin,p,cs.rowOwning,poses,cs.spanKind);
        tendon=min(tendon,r.design.tendon-span.requiredSpan_m);
    end
    if isAnkle && isfield(cs,'distalSleeveGuideRow')
        path=sum(vecnorm(diff(p(1:cs.distalSleeveGuideRow,:)),2,2));
        bay=min(bay,path-2*cs.fittingLength-max(m.sleeveLength_m));
    end
    for k=1:size(p,1)
        T=fr_chainT(kin,qHome,cs.rowOwning{k},'pelvis');
        z=T(3,1:3)*p(k,:)'+T(3,4);
        home=min(home,z);
    end
    lines{j}=struct('bodies',{cs.rowOwning},'points_m',p,'channel',r.channel, ...
        'count',cs.seed.nBPA,'gain',max(cs.txConfig.gain),'rest_mm',1000*r.design.L0, ...
        'tendon_mm',1000*r.design.tendon);
    for i=1:N
        q=poses{i};
        wire=worldPoints(kin,q,cs.rowOwning,p);
        x=primitives(C,cs,vessel{i},wire,j);
        parts{i}=[parts{i};x];
        if isAnkle
            start=cs.shankStopRow;
            foot=wire(start:min(start+1,size(wire,1)),:);
            minWire=min(minWire,segmentsAgainstBone(C,kin,q,[foot(1,:),foot(end,:),.002,2,0],{'tibia_r'}));
            T=fr_bodyT(kin,q,'tibia_r');
            solid=x(x(:,8)<2,:);
            a=(solid(:,1:3)-T(1:3,4)')*T(1:3,2);
            z=(solid(:,4:6)-T(1:3,4)')*T(1:3,2);
            footBay=min(footBay,min(min(a,z)-solid(:,7))+.425);
        end
    end
end
for i=1:N
    p=parts{i}; physical=p(p(:,8)<2,:);
    minBone=min(minBone,segmentsAgainstBone(C,kin,poses{i},physical,{'femur_r','tibia_r','talus_r'}));
end
signed=primary.*sign(targetPrimary);
capacity=sum(min(abs(targetPrimary),max(0,signed)))/max(sum(abs(targetPrimary)),eps);
D=struct('group',g,'lines',{lines},'parts',{parts},'bone_mm',1000*(minBone-C.clearanceBone_m), ...
    'foot_lead_mm',1000*minWire,'home_side_mm',1000*home,'peak_N',peak, ...
    'sleeve_bay_mm',1000*bay,'tendon_span_mm',1000*tendon,'foot_occupancy_mm',1000*footBay, ...
    'active',100*mean(hasForce),'common_active',100*mean(active), ...
    'capacity',100*capacity,'target_met',100*mean(hasForce & signed>=abs(targetPrimary)-1e-6), ...
    'primary_Nm',primary,'human_Nm',targetPrimary);
end

function [m,meta,cs]=atPoses(C,kin,r,branch,poses,isAnkle,target)
cs=r.caseDef;
if isAnkle
    G=fr_ankleCoupledSetup(C,kin,cs,2);
    cs=G.cs;
elseif ~isempty(target)
    G=fr_proximalCoupled(C,kin,cs,target);
    if ~isempty(G); cs=G.cs; end
end
cs.passes{1}.qStructs=poses;
qg=resolved(kin,r,cs,poses);
qg.pairedPoseGrid=true;
[m,meta]=fr_mechBuild(C,cs,r.design,kin,branch,qg);
end

function [t,a]=primaryCurve(C,kin,r,cs,branch,isAnkle)
coord=cs.passes{1}.sweptCoord;
poses=cs.passes{1}.qStructs;
cross=cs.primaryCrossing;
if isAnkle
    G=fr_ankleCoupledSetup(C,kin,r.caseDef,2);
    cs=G.cs;
    cross=1+strcmp(coord,G.second);
end
cs.passes{1}.sweptCoord=coord;
cs.passes{1}.qStructs=poses;
qg=resolved(kin,r,cs,poses);
qg.pairedPoseGrid=true;
m=fr_mechBuild(C,cs,r.design,kin,branch,qg);
a=transmission(m);
t=fr_coordinateTorque(cs,kin,m,1,cross);
t(~a|~isfinite(t))=0;
end

function qg=resolved(kin,r,cs,poses)
mask=cs.rowUse{1};
from=r.caseDef.rowOwning(mask); to=cs.rowFrame{1}(mask);
p=r.p(mask,:);
for j=1:numel(to)
    if strcmp(from{j},to{j}); continue; end
    T=fr_chainT(kin,kin.qDefault,from{j},to{j});
    p(j,:)=(T(1:3,1:3)*p(j,:)'+T(1:3,4))';
end
qg=struct('passIdx',1,'qList',{poses},'resolvedRoute',struct('points',p,'frames',{to}));
end

function a=transmission(m)
a=isfinite(m.Ftendon(:)) & m.Ftendon(:)>1e-6 & ...
    ~m.PulleySlack(:) & ~m.PulleyInfeasible(:) & ...
    m.strain_p(:)<=1-m.kmax/m.rest+1e-9;
end

function x=primitives(C,cs,v,wire,channel)
x=zeros(0,9);
for k=1:numel(v)
    x(end+1,:)=[v(k).p1,v(k).p2,v(k).r,0,100*channel+v(k).branch]; %#ok<AGROW>
end
[p,r,b]=fr_fittingSpheres(v,cs.fittingLength);
x=[x;p,p,r,ones(size(r)),100*channel+b];
for k=1:cs.seed.nBPA
    line=v([v.branch]==k);
    active=line([line.section]==1);
    tip=active(end).p2;
    base=v([v.branch]==1 & [v.section]==1);
    offset=tip-base(end).p2;
    id=100*channel+k;
    if cs.hasTackle
        family=C.families(cs.txConfig.tackleFamily);
        u=tip-active(end).p1; u=u/max(norm(u),eps);
        radius=(family.od_mm/2+2)/1000;
        x(end+1,:)=[tip,tip+.012*u,radius,1,id]; %#ok<AGROW>
    end
    wireStart=min(size(wire,1),max(2,cs.cross(1)));
    if isfield(cs,'shankStopRow')
        wireStart=cs.shankStopRow;
    end
    if any([line.section]==2)
        distal=line([line.section]==2);
        route=[tip;wire(wireStart:end-1,:)+offset;distal(end).p2];
    else
        route=[tip;wire(wireStart:end,:)+offset];
    end
    for j=1:size(route,1)-1
        x(end+1,:)=[route(j,:),route(j+1,:),.001,2,id]; %#ok<AGROW>
    end
    if isfield(cs,'shankStopRow') || strcmp(cs.spanKind,'knee')
        stop=wire(wireStart,:)+offset;
        x(end+1,:)=[stop,stop,.009,1,id]; %#ok<AGROW>
    end
end
end

function d=segmentsAgainstBone(C,kin,q,parts,bodies)
d=inf;
for b=bodies
    T=fr_bodyT(kin,q,b{1});
    for j=1:size(parts,1)
        p=parts(j,:);
        a=(T(1:3,1:3)'*(p(1:3)'-T(1:3,4)))';
        z=(T(1:3,1:3)'*(p(4:6)'-T(1:3,4)))';
        n=max(2,ceil(norm(z-a)/.005)+1);
        v=a+linspace(0,1,n)'*(z-a);
        delta=fr_surfaceDistance(C,b{1},v);
        d=min(d,min(delta)-norm(z-a)/(2*(n-1))-p(7));
    end
end
end

function p=worldPoints(kin,q,bodies,p)
for j=1:size(p,1)
    T=fr_bodyT(kin,q,bodies{j});
    p(j,:)=(T(1:3,1:3)*p(j,:)'+T(1:3,4))';
end
end

function L=routeLength(kin,p,b)
L=sum(vecnorm(diff(worldPoints(kin,kin.qDefault,b,p)),2,2));
end

function v=radical(i,p)
v=0; f=1/p;
while i>0
    v=v+f*mod(i,p); i=floor(i/p); f=f/p;
end
end

function write(file,value)
temp=[tempname(fileparts(file)) '.json'];
fid=fopen(temp,'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(value)); fclose(fid); movefile(temp,file,'f');
end
