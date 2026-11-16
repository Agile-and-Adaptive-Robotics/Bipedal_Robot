function fr_writeAnkleAssembly(C,kin,old,out,routes,solution,reference,selected)
% Publish model-frame paths and dense solved hardware for defense visuals.
T=solution.stats;
B=reference.stats;
for key={'ankle_met_pct','second_met_pct','ankle_coverage_pct','second_coverage_pct','no_force_pct'}
    field=key{1};
    T.(['before_' field])=B.(field);
end
writetable(T,fullfile(out,'torque_and_force.csv'));
writetable(solution.pairs,fullfile(out,'pair_clearance.csv'));
writetable(reference.pairs,fullfile(out,'before_pair_clearance.csv'));
actuators=cell(size(routes));
for j=1:numel(routes)
    r=routes{j}; d=solution.data{j};
    gain=max(r.caseDef.txConfig.gain);
    reeving=fr_travelReeving(gain);
    actuators{j}=struct('group',r.group,'channel',r.channel,'diameter_mm',20, ...
        'rest_mm',1000*r.design.L0,'tendon_mm',1000*r.design.tendon, ...
        'gain',gain,'sheaves',reeving.sheaveCount,'guideAdded',r.guideAdded, ...
        'offset_mm',1000*r.shift_m,'peakPerTendon_N',max(d.load,[],'omitnan'), ...
        'boneMargin_mm',1000*d.boneMargin_m,'leadMargin_mm',1000*d.leadMargin_m, ...
        'minimumHomeZ_mm',1000*fr_ankleHomeSide(kin,r), ...
        'failedConstraints',{d.fail},'bodies',{r.caseDef.rowOwning},'points_m',r.p);
end
S=struct('model',C.osimRobot,'anatomy',C.osimStock,'actuators',{actuators}, ...
    'selectedLayout',selected,'rigidConflicts',solution.rigidConflicts, ...
    'beforeRigidConflicts',reference.rigidConflicts,'localFailures',solution.localFailures, ...
    'beforeLocalFailures',reference.localFailures,'coverage',solution.coverage, ...
    'beforeCoverage',reference.coverage,'poseCount',size(solution.keys,1), ...
    'activeForce',1-mean(solution.stats.no_force_pct)/100, ...
    'beforeActiveForce',1-mean(reference.stats.no_force_pct)/100, ...
    'maxVesselOverlap_mm',max(solution.pairs.sleeve_overlap_mm),'old',old);
write(fullfile(out,'assembly_routes.json'),S);
source=fullfile(out,'joint_search.mat');
if ~isfile(source); source=fullfile(fileparts(out),'joint_search.mat'); end
W=load(source,'slots');
before=num2cell(W.slots);
for phase={'home','clash_pose'}
    if strcmp(phase{1},'home'); angle=[0 0 0]; else; angle=[20 35 20]; end
    [~,i]=min(vecnorm(solution.keys-angle,2,2));
    q=solution.data{1}.poses{i};
    bodies={'tibia_r','talus_r','calcn_r','toes_r'};
    transforms=struct();
    for j=1:numel(bodies)
        transforms.(bodies{j})=fr_chainT(kin,q,bodies{j},'tibia_r');
    end
    cloud=load(fullfile(C.boneCloudDir,'Tibia_Mesh_Points.txt'))*1000;
    scene=struct('angles_deg',solution.keys(i,:),'transforms',transforms, ...
        'anatomy',C.osimStock,'after',{sceneParts(kin,q,routes,solution,i,cloud)}, ...
        'before',{sceneParts(kin,q,before,reference,i,cloud)});
    write(fullfile(out,['scene_' phase{1} '.json']),scene);
end
save(fullfile(out,'signed_torque_grids.mat'),'solution','reference','-v7');
end

function parts=sceneParts(kin,q,routes,S,i,cloud)
parts=cell(size(routes));
T=fr_bodyT(kin,q,'tibia_r');
for j=1:numel(routes)
    d=S.data{j}; v=d.sleeves{i};
    capsules=cell(1,numel(v));
    for k=1:numel(v)
        capsules{k}=struct('points_mm',local([v(k).p1;v(k).p2],T),'radius_mm',1000*v(k).r);
    end
    [p,r]=fr_fittingSpheres(v,d.fit);
    extras=d.extra{i};
    p=[p;extras(:,1:3)];
    r=[r;extras(:,4)];
    first=v(1); last=v(end);
    u=(first.p2-first.p1)/max(eps,norm(first.p2-first.p1));
    at=local(first.p1-d.fit*u,T);
    [~,nearest]=min(vecnorm(cloud-at,2,2));
    if isfield(d,'paths')
        lead=[last.p2;d.paths{i}];
    else
        c=routes{j}.caseDef;
        p0=routes{j}.p;
        lead=last.p2;
        for k=c.distalSleeveGuideRow:size(p0,1)
            W=fr_bodyT(kin,q,c.rowOwning{k});
            lead(end+1,:)=(W(1:3,1:3)*p0(k,:)'+W(1:3,4))'; %#ok<AGROW>
        end
    end
    parts{j}=struct('group',routes{j}.group,'channel',routes{j}.channel, ...
        'vessels',{capsules},'hardware_mm',[local(p,T),1000*r], ...
        'wire_mm',local(lead,T),'mount_reach_mm',[cloud(nearest,:);at]);
end
end

function p=local(p,T)
p=1000*(T(1:3,1:3)'*(p'-T(1:3,4)))';
end

function write(file,value)
fid=fopen(file,'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(value)); fclose(fid);
end
