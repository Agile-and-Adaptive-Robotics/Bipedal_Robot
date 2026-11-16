function fr_exportGracilisPilotSweeps(outputDir)
% Evaluate the saved candidate at matching isolated knee and hip sweeps.
C = fr_setup();
kin = fr_kinematics(C.osimStock);
S = load(fullfile(outputDir,'grac_compliant.mat'),'result');
assert(~isempty(S.result.best),'No compliant gracilis length fit to plot');
c = S.result.best;
cs = c.caseDef;
mapFile = fullfile(C.repo,'Code','Matlab','Mesh_Optimization', ...
    'Routing27_v2','actuator_map_27.json');
map = jsondecode(fileread(mapFile));
coords = {'knee_angle_r','hip_adduction_r'};
rows = cell(1,2);
for p = 1:2
    coord = coords{p};
    if p==1
        limits = [-120 0];
    else
        limits = map.rom_table_deg.(coord);
    end
    angles = linspace(limits(1),limits(2),31);
    pass = cs.passes{p};
    pass.sweptCoord = coord;
    pass.th1 = deg2rad(angles);
    pass.qStructs = cell(size(angles));
    for j = 1:numel(angles)
        q = kin.qDefault;
        q.(coord) = pass.th1(j);
        pass.qStructs{j} = q;
    end
    cs.passes{p} = pass;
    qg = struct('passIdx',p,'qList',{pass.qStructs});
    m = fr_mechBuild(C,cs,c.design,kin,'compliant',qg);
    tau = fr_coordinateTorque(cs,kin,m,p,1);
    invalid = m.PulleyInfeasible | ~isfinite(tau);
    slack = logical(m.PulleySlack) & ~invalid;
    tau(invalid) = NaN;
    tau(slack) = 0;
    rows{p} = struct('coordinate',coord,'angles_deg',angles, ...
        'robot_Nm',tau,'slack',slack,'invalid',invalid,'poses',{pass.qStructs});
end
grid = coupledGrid(C,kin,cs,c.design,map.rom_table_deg.hip_adduction_r);
out = struct('humanModel',C.osimStock,'muscle','grac_r','branch','compliant', ...
    'sweeps',{rows},'coupled',grid);
fid = fopen(fullfile(outputDir,'gracilis_two_dof_sweeps.json'),'w');
assert(fid>=0);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(out));
end

function grid = coupledGrid(C,kin,cs,design,hipLimits)
% Native BiPulley Cartesian grid: T1 = hip, T2 = knee.
p = 2;
pass = cs.passes{p};
pass.T2 = @(k,q) fr_chainT(k,q,'tibia_r','femur_r');
bodies = cs.rowFrame{p}(cs.rowUse{p});
cs.crossFor{p} = [find(~strcmp(bodies,'pelvis'),1),find(strcmp(bodies,'tibia_r'),1)];
assert(isequal(unique(bodies,'stable'),{'pelvis','femur_r','tibia_r'}) || ...
       isequal(unique(bodies,'stable'),{'pelvis','tibia_r'}));
n = 11;
hip = linspace(hipLimits(1),hipLimits(2),n);
knee = linspace(-120,0,n);
qh = cell(1,n); qk = cell(1,n); poses = cell(1,n*n);
for i = 1:n
    qh{i} = kin.qDefault; qh{i}.hip_adduction_r = deg2rad(hip(i));
    qk{i} = kin.qDefault; qk{i}.knee_angle_r = deg2rad(knee(i));
end
for j = 1:n
    for i = 1:n
        q = qh{i}; q.knee_angle_r = qk{j}.knee_angle_r;
        poses{i+(j-1)*n} = q;
    end
end
cs.passes{p} = pass;
m = fr_mechBuild(C,cs,design,kin,'compliant',struct('passIdx',p,'qList',{qh},'qList2',{qk}));
[I,J] = ndgrid(1:n,1:n);
indices = [I(:),J(:)];
pass.qStructs = poses;
pass.sweptCoord = 'hip_adduction_r'; cs.passes{p} = pass;
hipTorque = fr_coordinateTorque(cs,kin,m,p,1,indices);
pass.sweptCoord = 'knee_angle_r'; cs.passes{p} = pass;
kneeTorque = fr_coordinateTorque(cs,kin,m,p,2,indices);
slack = logical(m.cell.PulleySlack(:));
invalid = logical(m.cell.PulleyInfeasible(:)) | ~isfinite(hipTorque) | ~isfinite(kneeTorque);
hipTorque(slack & ~invalid) = 0; kneeTorque(slack & ~invalid) = 0;
hipTorque(invalid) = NaN; kneeTorque(invalid) = NaN;
grid = struct('hip_deg',hip,'knee_deg',knee,'hip_Nm',hipTorque,'knee_Nm',kneeTorque, ...
    'slack',slack,'invalid',invalid,'poses',{poses});
end
