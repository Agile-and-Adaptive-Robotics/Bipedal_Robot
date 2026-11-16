function fr_fixAnkleHomeSide(old,first,out)
% Repack mounts that cross the home pelvis centerline on the right leg.
C=fr_setup(); kin=fr_kinematics(C.osimStock);
if ~isfolder(out); mkdir(out); end
S=load(fullfile(first,'assembled_validation.mat'),'routes','solution','reference');
routes=S.routes; groups=S.solution.stats.group';
mask=cellfun(@(r)fr_ankleHomeSide(kin,r)<0,routes);
n=3;
cache=fr_checkAssembledAnkleLayout(C,kin,routes,groups,n);
pool=[];
try
    cluster=parcluster('local');
    cluster.JobStorageLocation=fullfile(out,'parallel_jobs');
    if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
    pool=parpool(cluster,min(7,cluster.NumWorkers));
catch ME
    fprintf('Worker startup failed: %s. Continuing serially.\n',ME.message);
end
for index=find(mask)
    file=fullfile(out,sprintf('home_route_%02d.mat',index));
    if isfile(file)
        T=load(file,'routes','cache'); routes=T.routes; cache=T.cache;
        continue
    end
    [X,Y,Z]=ndgrid([-65,-45,-25,0,25,45,65],[-38,-18,0,20,38],[0,10]);
    moves=[X(:),Y(:),Z(:)]/1000;
    attempts=cell(1,size(moves,1));
    if isempty(pool)
        for j=1:numel(attempts)
            attempts{j}=probe(C,kin,routes,groups,n,cache,index,moves(j,:));
        end
    else
        parfor j=1:numel(attempts)
            attempts{j}=probe(C,kin,routes,groups,n,cache,index,moves(j,:));
        end
    end
    [value,j]=min(cellfun(@(r)r.score,attempts));
    assert(isfinite(value),'No eligible right-side mount placement');
    q=attempts{j}; routes=q.routes; cache=q.cache;
    save(file,'routes','cache','-v7');
    fprintf('%s %d home minimum +z %.1f mm. Assembly has %d coarse rigid conflicts.\n', ...
        routes{index}.group,routes{index}.channel,1000*fr_ankleHomeSide(kin,routes{index}), ...
        cache.rigidConflicts);
end
assert(all(cellfun(@(r)fr_ankleHomeSide(kin,r)>=0,routes)));
solution=fr_checkAssembledAnkleLayout(C,kin,routes,groups,7);
reference=S.reference; selected=1;
save(fullfile(out,'assembled_validation.mat'),'routes','solution','reference','selected','-v7');
fr_writeAnkleAssembly(C,kin,old,out,routes,solution,reference,selected);
fprintf('Dense right-side layout: %d conflicts, %d local failures, %.1f%% torque, %.1f%% force coverage.\n', ...
    solution.rigidConflicts,solution.localFailures,100*solution.coverage, ...
    100-mean(solution.stats.no_force_pct));
end

function result=probe(C,kin,routes,groups,n,cache,index,move)
r=routes{index}; cs=r.caseDef;
home=fr_ankleHomeSide(kin,r);
T=fr_chainT(kin,kin.qDefault,'tibia_r','pelvis');
shift=(move+[0,0,-home+.003])*T(1:3,1:3);
guide=cs.distalSleeveGuideRow;
stop=cs.shankStopRow;
r.p(1:guide,:)=r.p(1:guide,:)+shift;
r.p(guide+1:stop,:)=r.p(guide+1:stop,:)+.32*shift;
for k=1:stop
    z=T(3,1:3)*r.p(k,:)'+T(3,4);
    r.p(k,:)=r.p(k,:)+max(0,.003-z)*T(3,1:3);
end
r.shift_m=r.shift_m+shift;
if fr_ankleHomeSide(kin,r)<0
    result=struct('score',inf);
    return
end
routes{index}=r;
s=fr_checkAssembledAnkleLayout(C,kin,routes,groups,n,cache,index);
best=score(s);
if ~isempty(s.data{index}.fail)
    G=fr_resolveAnkleNative(C,kin,cs,r.p,4);
    r.design=fr_refitMountLengths(C,kin,cs,r.design,'compliant',G);
    Q=routes; Q{index}=r;
    v=fr_checkAssembledAnkleLayout(C,kin,Q,groups,n,cache,index);
    if score(v)<best; s=v; best=score(v); routes=Q; end
end
result=struct('score',best+.02*norm(move)*1000,'routes',{routes},'cache',s);
end

function s=score(a)
s=a.score+120*(1-a.coverage)+20*max(a.pairs.sleeve_overlap_mm)+5*a.localFailures;
end
