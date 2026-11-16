function fr_clearAnkleInterfaces(old,first,out,stagger)
% Move individual supports against the exact simultaneous assembly objective.
C=fr_setup(); kin=fr_kinematics(C.osimStock);
if nargin<4; stagger=false; end
if ~isfolder(out); mkdir(out); end
S=load(fullfile(first,'assembled_validation.mat'),'routes','solution','reference');
routes=S.routes; groups=S.solution.stats.group';
n=3;
cache=fr_checkAssembledAnkleLayout(C,kin,routes,groups,n);
pool=[];
try
    cluster=parcluster('local');
    cluster.JobStorageLocation=fullfile(out,'parallel_jobs');
    if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
    pool=parpool(cluster,min(8,cluster.NumWorkers));
catch ME
    fprintf('Worker startup failed: %s. Continuing serially.\n',ME.message);
end
for sweep=1:2
    f=fullfile(out,sprintf('sweep_%d.mat',sweep));
    if isfile(f)
        Q=load(f,'routes','cache'); routes=Q.routes; cache=Q.cache;
        continue
    end
    W=sortrows(cache.pairs,'rigid_margin_mm');
    W=W(~W.rigid_pass,:);
    for k=1:min(height(W),10)
        a=find(cellfun(@(r)strcmp(r.group,W.group_a{k}) && r.channel==W.channel_a(k),routes),1);
        b=find(cellfun(@(r)strcmp(r.group,W.group_b{k}) && r.channel==W.channel_b(k),routes),1);
        attempts=cell(1,16+16*stagger);
        if isempty(pool)
            for i=1:numel(attempts)
                attempts{i}=probe(C,kin,routes,groups,n,cache,[a b],i,stagger);
            end
        else
            parfor i=1:numel(attempts)
                attempts{i}=probe(C,kin,routes,groups,n,cache,[a b],i,stagger);
            end
        end
        best=cost(cache); winner=0;
        for i=1:numel(attempts)
            if attempts{i}.score<best
                best=attempts{i}.score; winner=i;
            end
        end
        if winner>0
            q=attempts{winner}; routes=q.routes; cache=q.cache;
        end
    end
    save(f,'routes','cache','-v7');
    fprintf('Exact interface sweep %d: %d conflicts, %d local failures, %.1f%% torque coverage.\n', ...
        sweep,cache.rigidConflicts,cache.localFailures,100*cache.coverage);
end
solution=fr_checkAssembledAnkleLayout(C,kin,routes,groups,7);
reference=S.reference;
if cost(solution)>=cost(S.solution)
    solution=S.solution;
    routes=S.routes;
end
selected=1;
save(fullfile(out,'assembled_validation.mat'),'routes','solution','reference','selected','-v7');
fr_writeAnkleAssembly(C,kin,old,out,routes,solution,reference,selected);
fprintf('Dense exact layout: %d rigid conflicts, %d local failures, %.1f%% torque coverage.\n', ...
    solution.rigidConflicts,solution.localFailures,100*solution.coverage);
end

function s=cost(a)
s=a.score+4*a.localFailures+150*(1-a.coverage)+4*max(a.pairs.sleeve_overlap_mm);
end

function result=probe(C,kin,routes,groups,n,cache,ids,j,stagger)
width=8+8*stagger;
which=1+floor((j-1)/width);
index=ids(which);
direction=[1,0,0;-1,0,0;0,0,1;0,0,-1;1,0,1;1,0,-1;-1,0,1;-1,0,-1];
r=routes{index};
cs=r.caseDef;
stop=cs.shankStopRow;
action=mod(j-1,width)+1;
if stagger && action<=8
    shift=.014*direction(action,:)/norm(direction(action,:));
    r.p(stop,:)=r.p(stop,:)+shift;
    if r.guideAdded; r.p(stop-1,:)=r.p(stop-1,:)+.5*shift; end
else
    if stagger
        axial=[0,.016,0;0,-.016,0;0,.028,0;0,-.028,0; ...
            .012,.018,0;-.012,.018,0;0,.018,.012;0,.018,-.012];
        step=axial(action-8,:);
    else
        step=.011*direction(action,:)/norm(direction(action,:));
    end
    r.p(1:cs.distalSleeveGuideRow,:)=r.p(1:cs.distalSleeveGuideRow,:)+step;
    r.p(cs.distalSleeveGuideRow+1:stop,:)=r.p(cs.distalSleeveGuideRow+1:stop,:)+.5*step;
    r.shift_m=r.shift_m+step;
end
routes{index}=r;
s=fr_checkAssembledAnkleLayout(C,kin,routes,groups,n,cache,index);
best=cost(s);
if ~strcmp(r.group,'ext_dig') && ~isempty(s.data{index}.fail)
    G=fr_resolveAnkleNative(C,kin,cs,r.p,4);
    r.design=fr_refitMountLengths(C,kin,cs,r.design,'compliant',G);
    rr=routes;
    rr{index}=r;
    check=fr_checkAssembledAnkleLayout(C,kin,rr,groups,n,cache,index);
    value=cost(check);
    if value<best
        s=check; best=value; routes=rr;
    end
end
result=struct('score',best,'routes',{routes},'cache',s);
end
