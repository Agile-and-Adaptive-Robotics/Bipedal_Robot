function fr_exportAnkleClearanceWitnesses(runDir)
% Reproduce the governing clearance samples from the optimized checkpoint.
C=fr_setup(); C.exactSegmentDistance=true; C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
names={'ext_dig_compliant','per_tert_compliant','tib_post_rigid','flex_dig_compliant'};
out=cell(1,numel(names));
for i=1:numel(names)
    name=names{i};
    S=load(fullfile(runDir,[name '.mat']),'result','record');
    c=S.result.chosen; branch=S.record.branch;
    [E,D]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,branch,11);
    b=D.bone; lead=D.lead;
    assert(abs(b.margin_m-S.record.coupled.boneMargin_m)<1e-9);
    assert(abs(lead.margin_m-S.record.coupled.leadMargin_m)<1e-9);
    b.ankle_deg=E.ankle_deg(b.poseIndex);
    b.second_deg=E.second_deg(b.poseIndex);
    cap=D.sleeves{b.poseIndex}(b.capsuleIndex);
    T=fr_bodyT(kin,E.poses{b.poseIndex},b.body);
    b.sleeve_m=(T(1:3,1:3)'*([cap.p1;cap.p2]'-T(1:3,4)))';
    dir=diff(b.sleeve_m); dir=dir/norm(dir);
    s=dot(b.closest_m-b.sleeve_m(1,:),dir);
    b.part='inflatable vessel';
    if s < 0; b.part='proximal fitting envelope'; end
    if s > norm(diff(b.sleeve_m)); b.part='distal fitting envelope'; end
    lead.ankle_deg=E.ankle_deg(lead.poseIndex);
    lead.second_deg=E.second_deg(lead.poseIndex);
    probes=localProbes(C,kin,c,branch,E);
    out{i}=struct('name',name,'bone',b,'lead',lead,'second',E.second, ...
        'diameter_mm',S.record.diameter_mm,'count',S.record.count,'gain',S.record.gain, ...
        'bundleMargin_m',E.bundleMargin_m,'tendonShortfall_m', ...
        max(0,E.tendonRequired_m-c.design.tendon),'probes',{probes});
    fprintf('%s: %s margin %.2f mm at %.1f / %.1f deg, lead %.2f mm at %.1f / %.1f deg\n', ...
        name,b.part,1000*b.margin_m,b.ankle_deg,b.second_deg, ...
        1000*lead.margin_m,lead.ankle_deg,lead.second_deg);
end
outDir=fullfile(runDir,'clearance_diagnostics');
if ~isfolder(outDir); mkdir(outDir); end
fid=fopen(fullfile(outDir,'CLEARANCE_WITNESSES.json'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(out,'PrettyPrint',true));
end

function out=localProbes(C,kin,c,branch,base)
% Bounded one-change probes; report new constraints instead of choosing
% a correction based on a single pose.
cs=c.caseDef; design=c.design;
out={};
origin=find([cs.movable.row]==1,1);
if ~isempty(origin)
    j=3*origin-2:3*origin;
    for coord=1:3
        for side=1:2
            bounds=[cs.lb(j(coord)),cs.ub(j(coord))];
            d=design;
            d.uv(origin,coord)=.5*(bounds(side)+design.uv(origin,coord));
            if d.uv(origin,coord)==design.uv(origin,coord); continue; end
            e=fr_evalCoupledAnkle(C,kin,cs,d,branch,11);
            out{end+1}=probe(sprintf('origin variable %d, side %d',coord,side),e,base,d); %#ok<AGROW>
        end
    end
end
foot=find(ismember(cs.rowOwning,{'talus_r','calcn_r','toes_r'}),1);
stop=foot-1;
if stop>1 && strcmp(cs.rowOwning{stop},'tibia_r') && ~any([cs.movable.row]==stop)
    patch=fr_patch(C,'Tibia','tibia_r',cs.seedLoc(stop,:),15);
    offsets=[0 0; .5 0; -.5 0; 0 .5; 0 -.5];
    for j=1:size(offsets,1)
        modified=cs;
        point=patch.project(offsets(j,1),offsets(j,2),C.standoff_mm(2));
        modified.seedLoc(stop,:)=point(:)';
        e=fr_evalCoupledAnkle(C,kin,modified,design,branch,11);
        out{end+1}=probe(sprintf('tibial stop patch offset %.1f %.1f',offsets(j,:)), ...
            e,base,design); %#ok<AGROW>
    end
end
end

function p=probe(label,e,base,d)
p=struct('change',label,'boneMargin_m',e.boneMargin_m,'leadMargin_m',e.leadMargin_m, ...
    'tendonShortfall_m',max(0,e.tendonRequired_m-d.tendon),'bundleMargin_m',e.bundleMargin_m, ...
    'badTransmission',nnz(e.slack|e.stroke|e.invalid), ...
    'maxAnkleTorqueChange_Nm',max(abs(e.ankle_Nm-base.ankle_Nm)), ...
    'maxSecondTorqueChange_Nm',max(abs(e.second_Nm-base.second_Nm)), ...
    'localPass',e.localPass);
end
