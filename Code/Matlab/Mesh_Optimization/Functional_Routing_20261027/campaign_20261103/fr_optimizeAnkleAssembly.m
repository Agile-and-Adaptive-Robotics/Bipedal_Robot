function fr_optimizeAnkleAssembly(old,out,workers,layoutCount)
% Coordinated packing of all ten ankle groups with independent load paths.
if nargin<3; workers=8; end
if nargin<4; layoutCount=8; end
C=fr_setup(); kin=fr_kinematics(C.osimStock);
C.exactSegmentDistance=true; C.fullFittingEnvelopes=true;
if ~isfolder(out); mkdir(out); end
groups={'soleus','tib_post','tib_ant','per_brev','per_long','per_tert', ...
    'flex_dig','flex_hal','ext_hal','ext_dig'};
grid=fullfile(out,'joint_search.mat');
if isfile(grid)
    S=load(grid,'slots','layouts','offsets','costs'); slots=S.slots;
    layouts=S.layouts; offsets=S.offsets; costs=S.costs;
else
    slots=separateRoutes(C,kin,old,groups);
    [layouts,offsets,costs]=coordinatedSearch(C,kin,slots,layoutCount);
    save(grid,'slots','layouts','offsets','costs','groups','-v7');
end
fprintf('Coordinated search: %d separate actuators, %d saved layout seeds.\n',numel(slots),size(layouts,1));
pool=[];
if workers>1
    try
        cluster=parcluster('local');
        cluster.JobStorageLocation=fullfile(out,'parallel_jobs');
        if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
        pool=parpool(cluster,min([workers,cluster.NumWorkers,numel(slots)]));
    catch ME
        fprintf('Worker startup failed: %s. Continuing serially.\n',ME.message);
    end
end
best=inf;
for k=1:size(layouts,1)
    file=fullfile(out,sprintf('layout_%02d.mat',k));
    if isfile(file)
        S=load(file,'solution'); solution=S.solution;
    else
        offsetsHere=offsets(layouts(k,:),:);
        routes=cell(1,numel(slots));
        if isempty(pool)
            for i=1:numel(slots)
                routes{i}=fitRoute(C,kin,slots(i),offsetsHere(i,:),out,k,i);
            end
        else
            parfor i=1:numel(slots)
                routes{i}=fitRoute(C,kin,slots(i),offsetsHere(i,:),out,k,i);
            end
        end
        solution=fr_checkAssembledAnkleLayout(C,kin,routes,groups,3);
        save(file,'solution','routes','-v7');
    end
    fprintf('Layout %d: %d rigid conflicts, %d local route failures, signed torque coverage %.1f%%.\n', ...
        k,solution.rigidConflicts,solution.localFailures,100*solution.coverage);
    if solution.score<best
        best=solution.score;
        selected=k;
    end
end
S=load(fullfile(out,sprintf('layout_%02d.mat',selected)),'routes');
routes=S.routes;
solution=fr_checkAssembledAnkleLayout(C,kin,routes,groups,7);
baseline=baselineRoutes(slots);
reference=fr_checkAssembledAnkleLayout(C,kin,baseline,groups,7);
save(fullfile(out,'assembled_validation.mat'),'routes','solution','reference','selected','-v7');
fr_writeAnkleAssembly(C,kin,old,out,routes,solution,reference,selected);
fprintf('Dense assembly: rigid conflicts %d -> %d, local route failures %d -> %d, coverage %.1f%% -> %.1f%%.\n', ...
    reference.rigidConflicts,solution.rigidConflicts,reference.localFailures,solution.localFailures, ...
    100*reference.coverage,100*solution.coverage);
end

function routes=baselineRoutes(slots)
routes=cell(1,numel(slots));
for i=1:numel(slots)
    s=slots(i);
    routes{i}=struct('caseDef',s.caseDef,'design',s.design,'p',s.p, ...
        'group',s.group,'channel',s.channel,'shift_m',[0 0 0],'guideAdded',false);
end
end

function slots=separateRoutes(C,kin,old,groups)
A=load(fullfile(old,'assembled_compliant.mat'),'caseData','keys');
P=load(fullfile(old,'edl_fixed_pair_separate_guides','independent_pair_dense.mat'),'c','record');
slots=struct('caseDef',{},'design',{},'p',{},'group',{},'channel',{},'proxy',{});
for j=1:numel(groups)
    name=groups{j};
    S=load(fullfile(old,[name '_compliant.mat']),'result','record');
    c=S.result.chosen;
    count=c.caseDef.seed.nBPA;
    if strcmp(name,'ext_dig'); c=P.c; count=2; end
    cs=c.caseDef;
    if ~isfield(cs,'shankStopRow')
        c=fr_prepareAnkleMountSearch(C,kin,c);
        cs=c.caseDef;
    end
    if strcmp(name,'ext_dig')
        p=P.record.routeParts{1}.points_m;
    else
        p=S.record.points_m;
        if size(p,1)~=numel(cs.rowOwning)
            G=fr_ankleCoupledSetup(C,kin,cs,2);
            [p,F]=fr_locForPass(C,G.cs,c.design,kin,1);
            for r=1:size(p,1)
                T=fr_chainT(kin,kin.qDefault,F{r},cs.rowOwning{r});
                p(r,:)=(T(1:3,1:3)*p(r,:)'+T(1:3,4))';
            end
        end
    end
    cs.seed.nBPA=1; cs.txConfig.nPulleyBPA=[1 1];
    cs.assemblySizing=true;
    a=A.caseData{j};
    [~,at]=min(vecnorm(A.keys,2,2));
    T=fr_bodyT(kin,a.poses{at},'tibia_r');
    sleeves=a.sleeves{at};
    entries=zeros(count,3);
    for b=1:count
        e=sleeves(find([sleeves.branch]==b,1));
        entries(b,:)=toTibia(e.p1,T);
    end
    shift=entries-mean(entries,1);
    for b=1:count
        pp=p;
        if strcmp(name,'ext_dig')
            pp=P.record.routeParts{b}.points_m;
        elseif count>1
            pp(1:cs.distalSleeveGuideRow,:)=pp(1:cs.distalSleeveGuideRow,:)+shift(b,:);
            pp(cs.shankStopRow,:)=pp(cs.shankStopRow,:)+.4*shift(b,:);
            for r=cs.shankStopRow+1:size(pp,1)
                Q=fr_chainT(kin,kin.qDefault,'tibia_r',cs.rowOwning{r});
                pp(r,:)=pp(r,:)+.18*shift(b,:)*Q(1:3,1:3)';
            end
        end
        proxy=proxyParts(C,a,b,cs,pp,T);
        slots(end+1)=struct('caseDef',cs,'design',c.design,'p',pp, ...
            'group',name,'channel',b,'proxy',proxy); %#ok<AGROW>
    end
end
end

function v=toTibia(v,T)
v=(T(1:3,1:3)'*(v'-T(1:3,4)))';
end

function P=proxyParts(C,a,channel,cs,p,T)
% Cover the longest solved vessel centerline with inflated sample spheres.
lengths=cellfun(@(v)sum(arrayfun(@(x)norm(x.p2-x.p1),v([v.branch]==channel))),a.sleeves);
[~,i]=max(lengths);
v=a.sleeves{i}; v=v([v.branch]==channel);
radius=max(cellfun(@(x)max([x([x.branch]==channel).r]),a.sleeves));
P=zeros(0,5);
for x=v
    n=max(2,ceil(norm(x.p2-x.p1)/.018)+1);
    ends=toTibia([x.p1;x.p2],T);
    points=ends(1,:)+linspace(0,1,n)'.*diff(ends);
    P=[P;points,repmat([radius,0],n,1)]; %#ok<AGROW>
end
[x,r]=fr_fittingSpheres(v,a.fit);
P=[P;toTibia(x,T),r,ones(size(r))];
if cs.hasTackle
    fam=C.families(cs.txConfig.tackleFamily);
    P=[P;toTibia(v(end).p2,T),(fam.od_mm/2+2)/1000,1];
end
P=[P;p(cs.shankStopRow,:),hypot(C.bowden.stopAllowance_mm,C.bowden.housingOD_mm)/2000,2];
end

function [layouts,offsets,values]=coordinatedSearch(C,kin,slots,layoutCount)
rng(8402,'twister');
[X,Z]=ndgrid(-.090:.015:.090);
offsets=[X(:),zeros(numel(X),1),Z(:)];
offsets=offsets(vecnorm(offsets,2,2)<=.108,:);
N=size(offsets,1); M=numel(slots);
unary=zeros(N,M);
centers=cell(1,M);
T=fr_chainT(kin,kin.qDefault,'tibia_r','pelvis');
shiftZ=offsets*T(3,1:3)';
for a=1:M
    p=slots(a).proxy;
    D=permute(offsets,[1 3 2]).*reshape(ones(size(p,1),1),1,[],1);
    D(:,p(:,5)==2,:)=.4*D(:,p(:,5)==2,:);
    D=D+permute(p(:,1:3),[3 1 2]);
    centers{a}=D;
    distance=reshape(fr_surfaceDistance(C,'tibia_r',reshape(D,[],3)),N,[]);
    deficit=max(0,p(:,4)'+C.clearanceBone_m-distance);
    stand=max(distance(:,p(:,5)~=2),[],2);
    unary(:,a)=2e4*max(deficit,[],2).^2 + .018*sum((offsets/.015).^2,2)+ ...
        5e3*max(0,stand-.100).^2;
    g=slots(a).caseDef.distalSleeveGuideRow;
    stop=slots(a).caseDef.shankStopRow;
    z=slots(a).p(1:stop,:)*T(3,1:3)'+T(3,4);
    outside=min(z(1:g))+shiftZ<.003 | min(z(g+1:stop))+.4*shiftZ<.003;
    unary(outside,a)=1e6;
end
pair=cell(M);
for a=1:M-1
    A=centers{a}; p=slots(a).proxy;
    for b=a+1:M
        B=centers{b}; q=slots(b).proxy;
        overlap=zeros(N); soft=zeros(N);
        same=strcmp(slots(a).group,slots(b).group);
        for i=1:size(p,1)
            for j=1:size(q,1)
                dx=A(:,i,1)-B(:,j,1)';
                dy=A(:,i,2)-B(:,j,2)';
                dz=A(:,i,3)-B(:,j,3)';
                hard=same || p(i,5)>0 || q(j,5)>0;
                gap=p(i,4)+q(j,4)+hard*C.clearanceEnvelope_m-sqrt(dx.^2+dy.^2+dz.^2);
                if hard; overlap=max(overlap,gap); else; soft=max(soft,gap); end
            end
        end
        pair{a,b}=1.1e4*overlap.^2+1.1*(overlap>.001)+1600*soft.^2;
        pair{b,a}=pair{a,b}';
    end
end
starts=180;
history=zeros(starts,M); scores=inf(starts,1);
for s=1:starts
    if s==1
        [~,x]=min(unary,[],1);
    else
        x=randi(N,1,M);
        if mod(s,3)==0
            [~,k]=min(scores(1:s-1));
            x=history(k,:);
            swap=randperm(M,ceil(M/3)); x(swap)=randi(N,1,numel(swap));
        end
    end
    for pass=1:28
        prior=x;
        for a=randperm(M)
            cost=unary(:,a);
            for b=[1:a-1 a+1:M]
                cost=cost+pair{a,b}(:,x(b));
            end
            [~,x(a)]=min(cost);
        end
        if isequal(prior,x); break; end
    end
    scores(s)=sum(unary(sub2ind(size(unary),x,1:M)));
    for a=1:M-1
        for b=a+1:M; scores(s)=scores(s)+pair{a,b}(x(a),x(b)); end
    end
    history(s,:)=x;
end
[~,order]=sort(scores);
[layouts,i]=unique(history(order,:),'rows','stable');
values=scores(order(i));
layouts=layouts(1:min(layoutCount,end),:); values=values(1:size(layouts,1));
fprintf('All-group packing seeds screened %d standoff sites and %d coordinated starts.\n',N,starts);
end

function r=fitRoute(C,kin,s,offset,out,layout,index)
file=fullfile(out,sprintf('layout_%02d_actuator_%02d.mat',layout,index));
if isfile(file); S=load(file,'r'); r=S.r; return; end
cs=s.caseDef;
base=s.p;
base(1:cs.distalSleeveGuideRow,:)=base(1:cs.distalSleeveGuideRow,:)+offset;
best=inf;
for fraction=[.25,.48,.72]
    p=base; p(cs.shankStopRow,:)=p(cs.shankStopRow,:)+fraction*offset;
    for via=0:1
        c=cs; path=p;
        if via
            k=c.shankStopRow-1;
            g=mean(path(k:k+1,:),1)+.15*offset;
            c=fr_addViaGuide(c,k,g);
            c.shankStopRow=c.shankStopRow+1;
            path=[path(1:k,:);g;path(k+1:end,:)];
        end
        for gain=unique([1,max(c.txConfig.gain)],'stable')
            probe=c;
            rr=fr_travelReeving(gain);
            probe.txConfig.gain=[gain,1];
            probe.txConfig.tackleLineParts=[rr.movingFalls,1];
            probe.hasTackle=gain>1;
            d=s.design;
            G=fr_resolveAnkleNative(C,kin,probe,path,4);
            if ~strcmp(s.group,'ext_dig')
                span=fr_spans(C,kin,G.qg.resolvedRoute.points,G.qg.resolvedRoute.frames, ...
                    G.cs.passes{1}.qStructs,'foot');
                d.tendon=max(d.tendon,span.requiredSpan_m+.001);
                d=fr_refitMountLengths(C,kin,probe,d,'compliant',G);
            end
            [e,details]=fr_evalCoupledAnkle(C,kin,probe,d,'compliant',4,G);
            local=fr_ankleLayoutLocal(probe,d,e);
            check=struct('p',path,'caseDef',probe);
            if fr_ankleHomeSide(kin,check)<0; local=local+1e5; end
            baselineGain=max(cs.txConfig.gain);
            score=local+2e-5*norm(path(probe.shankStopRow,:)-s.p(cs.shankStopRow,:))*1000 ...
                +.001*via+.002*(gain~=baselineGain);
            if score<best
                best=score;
                r=struct('caseDef',probe,'design',d,'p',path,'group',s.group, ...
                    'channel',s.channel,'shift_m',offset,'guideAdded',logical(via), ...
                    'sampledForce_N',mean(e.tendonForce_N,'omitnan'),'bone',details.bone);
            end
            if gain==1 && e.localPass; break; end
        end
    end
end
tmp=[tempname(out) '.mat']; save(tmp,'r','-v7'); movefile(tmp,file);
end
