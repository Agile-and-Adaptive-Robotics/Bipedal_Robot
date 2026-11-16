function fr_refineAnkleAssembly(old,first,out,workers)
% Close the local cable/length failures while preserving whole-layout fit.
if nargin<4; workers=8; end
C=fr_setup(); kin=fr_kinematics(C.osimStock);
if ~isfolder(out); mkdir(out); end
validated=fullfile(first,'assembled_validation.mat');
if isfile(validated)
    S=load(validated,'routes','solution','reference');
else
    layouts=dir(fullfile(first,'layout_??.mat'));
    score=inf;
    for i=1:numel(layouts)
        q=load(fullfile(first,layouts(i).name),'routes','solution');
        if q.solution.score<score; score=q.solution.score; S=q; end
    end
end
base=S.routes;
groups=S.solution.stats.group';
result=cell(size(base));
cluster=parcluster('local');
cluster.JobStorageLocation=fullfile(out,'parallel_jobs');
if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
pool=[];
try
    pool=parpool(cluster,min([workers,cluster.NumWorkers,numel(base)]));
catch ME
    fprintf('Worker startup failed: %s. Continuing serially.\n',ME.message);
end
if isempty(pool)
    for i=1:numel(base)
        result{i}=alternatives(C,kin,base{i},S.solution.data{i},out,i);
    end
else
    parfor i=1:numel(base)
        result{i}=alternatives(C,kin,base{i},S.solution.data{i},out,i);
    end
end

% Both passes score all actuators together after the local length fits.
layouts={base,cellfun(@(a)a{1},result,'UniformOutput',false)};
for round=1:4
    r=layouts{2};
    for i=1:numel(r)
        a=result{i};
        if numel(a)<2; continue; end
        if mod(i+round,4)==0; r{i}=a{2}; end
    end
    layouts{end+1}=r; %#ok<AGROW>
end
best=inf;
for j=1:numel(layouts)
    routes=layouts{j};
    candidate=fr_checkAssembledAnkleLayout(C,kin,routes,groups,3);
    score=candidate.score+5*candidate.localFailures;
    fprintf('Refinement %d: %d rigid conflicts, %d local failures, %.1f%% signed torque coverage.\n', ...
        j,candidate.rigidConflicts,candidate.localFailures,100*candidate.coverage);
    save(fullfile(out,sprintf('refinement_%02d.mat',j)),'routes','candidate','-v7');
    if score<best; best=score; selected=j; end
end
routes=layouts{selected};
solution=fr_checkAssembledAnkleLayout(C,kin,routes,groups,7);
baseline=load(validated,'routes','solution','reference');
reference=baseline.reference;
if solution.score+5*solution.localFailures>baseline.solution.score+5*baseline.solution.localFailures
    routes=baseline.routes;
    solution=baseline.solution;
    selected=1;
end
save(fullfile(out,'assembled_validation.mat'),'routes','solution','reference','selected','-v7');
fr_writeAnkleAssembly(C,kin,old,out,routes,solution,reference,selected);
fprintf('Dense refined assembly: rigid conflicts %d -> %d, local failures %d -> %d, coverage %.1f%% -> %.1f%%.\n', ...
    reference.rigidConflicts,solution.rigidConflicts,reference.localFailures,solution.localFailures, ...
    100*reference.coverage,100*solution.coverage);
end

function routes=alternatives(C,kin,seed,previous,out,index)
file=fullfile(out,sprintf('route_%02d.mat',index));
if isfile(file); S=load(file,'routes'); routes=S.routes; return; end
cs=seed.caseDef;
stop=cs.shankStopRow;
near=load(fullfile(C.boneCloudDir,'Tibia_Mesh_Points.txt'));
[~,j]=min(vecnorm(near-seed.p(stop,:),2,2));
outward=seed.p(stop,:)-near(j,:); outward(2)=0;
outward=outward/max(norm(outward),eps);
dest=fr_chainT(kin,kin.qDefault,cs.rowOwning{stop+1},'tibia_r');
foot=(dest(1:3,1:3)*seed.p(stop+1,:)'+dest(1:3,4))';
toward=foot-seed.p(stop,:);
toward=toward/max(norm(toward),eps);
moves=[0 0 0; .008*outward; .018*outward; .010*toward; .020*toward; ...
    -.010*toward; .014*outward+.012*toward; .015*outward-.012*toward; ...
    .012 0 0; -.012 0 0; 0 0 .012; 0 0 -.012; 0 .012 0; 0 -.012 0];
if isempty(previous.fail)
    moves=moves([1,2,4,6],:);
end
ankle=fr_chainT(kin,kin.qDefault,'talus_r','tibia_r');
list={}; values=[];
for j=1:size(moves,1)
    r=seed;
    r.p(stop,:)=r.p(stop,:)+moves(j,:);
    if r.p(stop,2)<ankle(2,4)+.003; continue; end
    candidates={r};
    if ~seed.guideAdded && j<=8 && ~isempty(previous.fail)
        g=stop-1;
        q=mean(r.p(g:g+1,:),1)+.014*outward;
        u=r;
        u.p=[r.p(1:g,:);q;r.p(g+1:end,:)];
        u.caseDef=fr_addViaGuide(cs,g,q);
        u.caseDef.shankStopRow=stop+1;
        u.guideAdded=true;
        candidates{end+1}=u; %#ok<AGROW>
    end
    for k=1:numel(candidates)
        r=candidates{k};
        G=fr_resolveAnkleNative(C,kin,r.caseDef,r.p,4);
        if ~strcmp(r.group,'ext_dig')
            variants=lengthProposals(C,kin,r,G);
        else
            variants=r.design;
        end
        for v=variants
            r.design=v;
            e=fr_evalCoupledAnkle(C,kin,r.caseDef,r.design,'compliant',4,G);
            [f,~]=fr_ankleLayoutLocal(r.caseDef,r.design,e);
            % Among locally comparable fits, preserve force and torque.
            force=e.tendonForce_N; force(e.noForce)=0;
            retained=max(0,mean(force,'omitnan'))/1000;
            score=f+.9*norm(moves(j,:))+.004*(r.guideAdded~=seed.guideAdded)-.05*retained;
            list{end+1}=r; values(end+1)=score; %#ok<AGROW>
        end
    end
end
[~,I]=sort(values);
routes=list(I(1:min(2,numel(I))));
if isempty(routes); routes={seed}; end
tmp=[tempname(out) '.mat']; save(tmp,'routes','-v7'); movefile(tmp,file);
end

function designs=lengthProposals(C,kin,r,G)
% Screen tendon and rest lengths against the actual equilibrium endpoints.
cs=r.caseDef; d=r.design;
m=fr_mechBuild(C,G.cs,d,kin,'compliant',G.qg); x=m.cell;
P=G.qg.resolvedRoute.points; F=G.qg.resolvedRoute.frames;
span=fr_spans(C,kin,P,F,G.cs.passes{1}.qStructs,'foot');
bay=sum(vecnorm(diff(P(1:cs.distalSleeveGuideRow,:)),2,2))-2*cs.fittingLength;
index=m.ActiveCrossing;
gain=x.PulleyGain(index);
lengths=x.MuscleLength(:);
delta=x.deltaL(:,:,index);
cap=1-C.LminFraction.d20;
minimum=range(lengths+(1-1/gain)*delta(:))/cap;
upper=max(.03,bay/(1-.08));
lower=max(.028,min(upper,minimum*.92));
rest=unique(max(lower,min(upper,[d.L0,linspace(lower,upper,7)])));
spring=x.kSpr*d.tendon;
deflect=x.DeflClosure_cAdd(:)-gain*x.cbBracket(:);
arc=0;
if isfield(cs,'bpaBendArc_m'); arc=C.Xi.Xi3*cs.bpaBendArc_m; end
rows=zeros(0,3);
for L=rest
    t0=max(span.requiredSpan_m,max(lengths)-2*cs.fittingLength-L-.045);
    t1=max(t0+.004,max(lengths)-2*cs.fittingLength-L+.085);
    t=unique([linspace(t0,t1,120),max(t0,d.tendon)]);
    initial=L+t+x.Xi0+2*cs.fittingLength-lengths;
    bend=arc*max(0,1-initial/(cap*L)).^2;
    required=gain*(initial+bend)-(gain-1)*delta(:);
    force=maxBPAforce(L,'20',C.pressure);
    coeff=deflect-t/(gain*spring);
    lo=-required+coeff*(force*festo4(20,0,C.pressure));
    hi=gain*cap*L-required+coeff*(force*festo4(20,1,C.pressure));
    wrong=mean(lo>0 | hi<0,1);
    reserve=mean(max(0,min(1,-lo./max(hi-lo,eps))),1);
    penalty=300*wrong+reserve+max(0,L-bay)*50;
    [value,k]=min(penalty);
    rows(end+1,:)=[value,L,t(k)]; %#ok<AGROW>
end
rows=sortrows(rows,1);
designs=repmat(d,1,min(2,size(rows,1))+1);
for i=2:numel(designs)
    designs(i).L0=rows(i-1,2);
    designs(i).tendon=rows(i-1,3);
end
end
