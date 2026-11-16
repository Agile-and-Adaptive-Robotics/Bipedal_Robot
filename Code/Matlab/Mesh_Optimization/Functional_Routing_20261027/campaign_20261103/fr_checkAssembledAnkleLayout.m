function S=fr_checkAssembledAnkleLayout(C,kin,routes,groups,n,previous,changed)
% All physical tendons and vessels evaluated at common ankle/foot poses.
rom=jsondecode(fileread(fullfile(C.repo,'Code','Matlab','Mesh_Optimization', ...
    'Routing27_v2','actuator_map_27.json'))).rom_table_deg;
ax=linspace(rom.ankle_angle_r(1),rom.ankle_angle_r(2),n);
sx=linspace(rom.subtalar_angle_r(1),rom.subtalar_angle_r(2),n);
mx=linspace(rom.mtp_angle_r(1),rom.mtp_angle_r(2),n);
[A,B,D]=ndgrid(ax,sx,mx);
keys=sortrows([A(:),B(:),D(:)]);
data=cell(size(routes)); local=0; localCost=0;
if nargin>=6
    data=previous.data;
else
    changed=1:numel(routes);
end
for i=1:numel(routes)
    r=routes{i};
    if ismember(i,changed); data{i}=evaluate(C,kin,r,keys,ax,sx,mx,n); end
    local=local+~isempty(data{i}.fail);
    localCost=localCost+data{i}.localCost;
end
pairs=cell(nchoosek(numel(routes),2),1);
k=0;
for i=1:numel(routes)-1
    for j=i+1:numel(routes)
        k=k+1;
        if nargin>=6 && ~ismember(i,changed) && ~ismember(j,changed)
            pairs{k}=table2struct(previous.pairs(k,:));
            continue
        end
        a=data{i}; b=data{j};
        assert(isequal(a.poses,b.poses));
        L=fr_compareAnkleLayers(C,kin,a.poses,a.sleeves,b.sleeves,a.fit,b.fit,false);
        hardware=min([L.rigidMargin_m,hardwareMargin(C,a,b),hardwareMargin(C,b,a)]);
        same=strcmp(routes{i}.group,routes{j}.group);
        if same
            hardware=min(hardware,-L.sleeveOverlap_m);
            hardware=min(hardware,wireWitness(a,b));
        end
        pairs{k}=struct('group_a',routes{i}.group,'channel_a',routes{i}.channel, ...
            'group_b',routes{j}.group,'channel_b',routes{j}.channel, ...
            'rigid_margin_mm',1000*hardware,'sleeve_overlap_mm',1000*L.sleeveOverlap_m, ...
            'rigid_pass',hardware>=-1e-4);
    end
end
T=struct2table(vertcat(pairs{:}));
[stats,coverage,torques]=scoreTorque(C,routes,groups,data,keys);
rigid=nnz(~T.rigid_pass);
depth=sum(max(0,-T.rigid_margin_mm));
soft=sum(T.sleeve_overlap_mm);
score=80*rigid+1.5*depth+2*localCost+50*(1-coverage)+.05*soft;
S=struct('keys',keys,'data',{data},'pairs',T,'stats',stats,'torques',{torques}, ...
    'coverage',coverage,'rigidConflicts',rigid,'localFailures',local,'score',score);
end

function s=evaluate(C,kin,r,keys,ax,sx,mx,n)
cs=r.caseDef;
isToe=any(strcmp(cs.rowOwning,'toes_r'));
if isToe
    fixed=sx; sweep=mx; other='subtalar_angle_r';
else
    fixed=mx; sweep=sx; other='mtp_angle_r';
end
count=size(keys,1);
pos=cell(1,count); caps=pos; extra=pos; wire=pos; paths=pos;
ankle=zeros(count,1); second=ankle; load=ankle; bad=false(count,1);
cells=zeros(size(keys)); fail={}; localCost=0; d0=[];
minBone=inf; minLead=inf;
for i=1:n
    K=kin; K.qDefault.(other)=deg2rad(fixed(i));
    G=fr_resolveAnkleNative(C,K,cs,r.p,n,{ax,sweep});
    [e,d]=fr_evalCoupledAnkle(C,K,cs,r.design,'compliant',n,G);
    [penalty,flags]=fr_ankleLayoutLocal(cs,r.design,e);
    minBone=min(minBone,e.boneMargin_m); minLead=min(minLead,e.leadMargin_m);
    fail=union(fail,flags); localCost=max(localCost,penalty);
    ids=(i-1)*n^2+(1:n^2);
    pos(ids)=e.poses; caps(ids)=d.sleeves;
    ankle(ids)=e.ankle_Nm; second(ids)=e.second_Nm; bad(ids)=e.noForce;
    load(ids)=e.tendonForce_N;
    cells(ids,1)=e.ankle_deg;
    if isToe
        cells(ids,2:3)=[repmat(fixed(i),numel(ids),1),e.second_deg(:)];
    else
        cells(ids,2:3)=[e.second_deg(:),repmat(fixed(i),numel(ids),1)];
    end
    for j=1:numel(ids)
        q=e.poses{j};
        p=toWorld(K,q,cs.rowOwning,r.p);
        wire{ids(j)}=p(cs.shankStopRow:end,:);
        paths{ids(j)}=p(cs.distalSleeveGuideRow:end,:);
        stopRadius=hypot(C.bowden.stopAllowance_mm,C.bowden.housingOD_mm)/2000;
        spheres=[p(cs.shankStopRow,:),stopRadius];
        if r.guideAdded
            spheres(end+1,:)=[p(cs.shankStopRow-1,:),.004]; %#ok<AGROW>
        end
        if cs.hasTackle
            family=C.families(cs.txConfig.tackleFamily);
            spheres(end+1,:)=[d.sleeves{j}(end).p2,(family.od_mm/2+2)/1000]; %#ok<AGROW>
        end
        extra{ids(j)}=spheres;
    end
    if isempty(d0) || abs(fixed(i))==min(abs(fixed)); d0=d; end
end
[values,order]=sortrows(cells);
assert(isequal(values,keys));
home=fr_ankleHomeSide(kin,r);
if home<0
    fail=union(fail,{'home attachment side'});
    localCost=localCost+100+1000*abs(home);
end
s=struct('poses',{pos(order)},'sleeves',{caps(order)},'extra',{extra(order)}, ...
    'wire',{wire(order)},'paths',{paths(order)},'detail',d0,'fit',cs.fittingLength,'ankle',ankle(order), ...
    'second',second(order),'load',load(order),'bad',bad(order),'isToe',isToe, ...
    'localCost',localCost,'fail',{fail},'boneMargin_m',minBone,'leadMargin_m',minLead, ...
    'homeZ_m',home);
end

function p=toWorld(kin,q,bodies,p)
for r=1:size(p,1)
    T=fr_bodyT(kin,q,bodies{r});
    p(r,:)=(T(1:3,1:3)*p(r,:)'+T(1:3,4))';
end
end

function margin=wireWitness(a,b)
margin=inf;
for i=1:numel(a.wire)
    A=a.wire{i}; B=b.wire{i};
    for j=1:size(A,1)-1
        for k=1:size(B,1)-1
            margin=min(margin,fr_segmentDistance(A(j,:),A(j+1,:),B(k,:),B(k+1,:))-.002);
        end
    end
end
end

function margin=hardwareMargin(C,a,b)
margin=inf;
for i=1:numel(a.poses)
    [p,r]=fr_fittingSpheres(b.sleeves{i},b.fit);
    targets=[b.extra{i};p,r];
    for h=a.extra{i}'
        q=h(1:3)';
        margin=min(margin,min(vecnorm(targets(:,1:3)-q,2,2)-targets(:,4)-h(4))-C.clearanceEnvelope_m);
        for v=b.sleeves{i}
            gap=fr_segmentDistance(q,q,v.p1,v.p2)-h(4)-v.r-C.clearanceEnvelope_m;
            margin=min(margin,gap);
        end
    end
end
end

function [T,total,torques]=scoreTorque(C,routes,groups,data,keys)
persistent source human
base=fullfile(C.camp,'run_20261106_ankle20_independent');
if isempty(source) || ~strcmp(source,base)
    human=cell(size(groups));
    for j=1:numel(groups)
        g=groups{j};
        if strcmp(g,'soleus')
            human{j}=jsondecode(fileread(fullfile(base,'soleus_compliant.json')));
        else
            human{j}=readtable(fullfile(base,[g '_compliant_coupled_torque.csv']));
        end
    end
    source=base;
end
rows=cell(size(groups)); torques=rows; allCoverage=[];
for j=1:numel(groups)
    g=groups{j};
    ids=find(cellfun(@(x)strcmp(x.group,g),routes));
    D=data(ids);
    first=zeros(size(keys,1),1); second=first; bad=true(size(first));
    for d=D
        x=d{1};
        a=x.ankle; b=x.second;
        a(x.bad)=0; b(x.bad)=0;
        first=first+a; second=second+b;
        bad=bad & x.bad;
    end
    if strcmp(g,'soleus')
        H=human{j};
        A=interp1(H.angles_deg,H.human_Nm,keys(:,1),'linear','extrap');
        B=nan(size(A)); hit2=NaN; cov2=NaN;
    else
        H=human{j};
        ax=unique(H.ankle_deg); sx=unique(H{:,2});
        F=griddedInterpolant({ax,sx},reshape(H.human_ankle_Nm,numel(ax),numel(sx)));
        cross=2+D{1}.isToe;
        A=F(keys(:,1),keys(:,cross));
        F.Values=reshape(H.human_second_Nm,numel(ax),numel(sx));
        B=F(keys(:,1),keys(:,cross));
        [hit2,cov2]=metrics(second,B,bad);
    end
    [hit,cov]=metrics(first,A,bad);
    allCoverage=[allCoverage,cov,cov2]; %#ok<AGROW>
    rows{j}=struct('group',g,'ankle_met_pct',100*hit,'second_met_pct',100*hit2, ...
        'ankle_coverage_pct',100*cov,'second_coverage_pct',100*cov2, ...
        'no_force_pct',100*mean(bad),'actuators',numel(ids), ...
        'local_failures',sum(cellfun(@(x)~isempty(x.fail),D)));
    torques{j}=struct('group',g,'ankle_Nm',first,'second_Nm',second, ...
        'humanAnkle_Nm',A,'humanSecond_Nm',B,'noForce',bad);
end
T=struct2table(vertcat(rows{:}));
total=mean(allCoverage,'omitnan');
end

function [hit,coverage]=metrics(t,h,bad)
p=t.*sign(h);
hit=mean(p>=abs(h)-1e-6 & ~bad);
coverage=sum(min(abs(h),max(0,p)))/max(sum(abs(h)),eps);
end
