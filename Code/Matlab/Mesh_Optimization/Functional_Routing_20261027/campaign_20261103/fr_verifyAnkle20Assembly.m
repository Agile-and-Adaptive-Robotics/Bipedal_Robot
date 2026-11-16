function fr_verifyAnkle20Assembly(folder,branches)
% Recheck all selected sleeves and fittings on simultaneous distal poses.
C=fr_setup();
C.exactSegmentDistance=true; C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
rom=jsondecode(fileread(fullfile(C.repo,'Code','Matlab','Mesh_Optimization', ...
    'Routing27_v2','actuator_map_27.json'))).rom_table_deg;
groups={'soleus','tib_post','tib_ant','per_brev','per_long','per_tert', ...
    'flex_dig','flex_hal','ext_hal','ext_dig'};
if nargin<2; branches={'rigid','compliant'}; end
for branch=branches
    name=branch{1};
    file=fullfile(folder,['assembled_' name '.mat']);
    if isfile(file); continue; end
    n=5;
    a=rom.ankle_angle_r; b=rom.subtalar_angle_r; d=rom.mtp_angle_r;
    [A,B,D]=ndgrid(linspace(a(1),a(2),n),linspace(b(1),b(2),n),linspace(d(1),d(2),n));
    [keys,~]=sortrows([A(:),B(:),D(:)]);
    caseData=cell(size(groups));
    for i=1:numel(groups)
        if strcmp(name,'compliant') && strcmp(groups{i},'ext_dig')
            S=load(fullfile(folder,'edl_fixed_pair_separate_guides', ...
                'independent_pair_dense.mat'),'c','record');
            x=collect(C,kin,S.c,name,keys,n,S.record.routeParts{1});
            y=collect(C,kin,S.c,name,keys,n,S.record.routeParts{2});
            for j=1:numel(x.poses)
                second=y.sleeves{j};
                for k=1:numel(second); second(k).branch=2; end
                x.sleeves{j}=[x.sleeves{j},second];
                x.hardware{j}=[x.hardware{j};y.hardware{j}];
            end
            caseData{i}=x;
        else
            S=load(fullfile(folder,[groups{i} '_' name '.mat']),'result');
            caseData{i}=collect(C,kin,S.result.chosen,name,keys,n);
        end
    end
    checks=cell(1,nchoosek(numel(groups),2));
    k=0;
    for a=1:numel(groups)-1
        for b=a+1:numel(groups)
            x=caseData{a}; y=caseData{b};
            assert(isequal(x.poses,y.poses));
            L=fr_compareAnkleLayers(C,kin,x.poses,x.sleeves,y.sleeves,x.fit,y.fit,false);
            rigid=min([L.rigidMargin_m,extrasMargin(C,x,y),extrasMargin(C,y,x)]);
            k=k+1;
            checks{k}=struct('group_a',groups{a},'group_b',groups{b}, ...
                'fitting_margin_mm',1000*L.rigidMargin_m, ...
                'rigid_margin_mm',1000*rigid, ...
                'vessel_overlap_mm',1000*L.sleeveOverlap_m, ...
                'fitting_pass',L.rigidMargin_m>=-1e-4,'rigid_pass',rigid>=-1e-4);
        end
    end
    T=struct2table(vertcat(checks{:}));
    writetable(T,fullfile(folder,['assembled_' name '.csv']));
    save(file,'T','caseData','keys','-v7');
    fprintf('%s assembled: %d/%d group pairs clear the modeled rigid hardware.\n', ...
        name,nnz(T.rigid_pass),height(T));
end
end

function d=collect(C,kin,c,branch,keys,n,nativeRoute)
cs=c.caseDef;
ankle=unique(keys(:,1));
subtalar=unique(keys(:,2));
mtp=unique(keys(:,3));
if any(strcmp(cs.rowOwning,'toes_r'))
    second=mtp; fixed=subtalar; other='subtalar_angle_r';
else
    second=subtalar; fixed=mtp; other='mtp_angle_r';
end
pose=cell(1,size(keys,1));
caps=pose;
hardware=pose;
list=zeros(size(keys));
for i=1:numel(fixed)
    K=kin;
    K.qDefault.(other)=deg2rad(fixed(i));
    G=fr_ankleCoupledSetup(C,K,cs,n,{ankle,second});
    if nargin>=7
        points=nativeRoute.points_m;
        frames=G.cs.rowFrame{1};
        for r=1:size(points,1)
            T=fr_chainT(K,K.qDefault,nativeRoute.bodies{r},frames{r});
            points(r,:)=(T(1:3,1:3)*points(r,:)'+T(1:3,4))';
        end
        G.qg.resolvedRoute=struct('points',points,'frames',{frames});
    end
    [e,details]=fr_evalCoupledAnkle(C,K,cs,c.design,branch,n,G);
    ids=(i-1)*n^2+(1:n^2);
    pose(ids)=e.poses;
    caps(ids)=details.sleeves;
    for j=1:numel(ids)
        spheres=zeros(0,4);
        if isfield(cs,'shankStopRow')
            r=cs.shankStopRow;
            T=fr_bodyT(K,e.poses{j},details.frames{r});
            point=(T(1:3,1:3)*details.loc(r,:)'+T(1:3,4))';
            radius=hypot(C.bowden.stopAllowance_mm,C.bowden.housingOD_mm)/2000;
            spheres(end+1,:)=[point,radius]; %#ok<AGROW>
        end
        if cs.hasTackle
            family=C.families(cs.txConfig.tackleFamily);
            for branchID=unique([details.sleeves{j}.branch])
                part=details.sleeves{j}([details.sleeves{j}.branch]==branchID);
                spheres(end+1,:)=[part(end).p2,(family.od_mm/2+2)/1000]; %#ok<AGROW>
            end
        end
        hardware{ids(j)}=spheres;
    end
    list(ids,1)=e.ankle_deg;
    if strcmp(other,'subtalar_angle_r')
        list(ids,2:3)=[repmat(fixed(i),numel(ids),1),e.second_deg(:)];
    else
        list(ids,2:3)=[e.second_deg(:),repmat(fixed(i),numel(ids),1)];
    end
end
[values,order]=sortrows(list);
assert(isequal(values,keys));
d=struct('poses',{pose(order)},'sleeves',{caps(order)}, ...
    'hardware',{hardware(order)},'fit',cs.fittingLength);
end

function margin=extrasMargin(C,a,b)
margin=inf;
for i=1:numel(a.poses)
    [centers,radii]=fr_fittingSpheres(b.sleeves{i},b.fit);
    targets=[b.hardware{i};centers,radii];
    for sphere=a.hardware{i}'
        p=sphere(1:3)';
        radius=sphere(4);
        if ~isempty(targets)
            gap=vecnorm(targets(:,1:3)-p,2,2)-targets(:,4)-radius;
            margin=min(margin,min(gap)-C.clearanceEnvelope_m);
        end
        for v=b.sleeves{i}
            gap=fr_segmentDistance(p,p,v.p1,v.p2)-radius-v.r;
            margin=min(margin,gap-C.clearanceEnvelope_m);
        end
    end
end
end
