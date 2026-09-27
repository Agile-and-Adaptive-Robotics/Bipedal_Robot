function sns_bone_pipeline_20260926()
%% sns_bone_pipeline_20260926  OpenSim bone geometry -> Simscape parts.
%% 1. curate the gait2392 bone STLs into Solid_Models\Simscape_Part_Library\OpenSim_Bones
%% 2. build a demo model: World -> spherical hip -> femur -> revolute knee
%%    -> tibia -> revolute ankle -> foot, all File Solid bones at Onyx density
%% 3. compile + 0.05 s sim (proves the route end to end)
%% NOTE: joint anchors are at body-frame origins for the demo; exact OpenSim
%% joint offsets can be pulled from the converted MJCF (gait2392_robot.mjc)
%% when the bones go into a real walker.

sns = fileparts(fileparts(mfilename('fullpath')));           % ...\Code\Matlab\SNS_Simscape
repo = fileparts(fileparts(fileparts(sns)));                 % ...\Bipedal_Robot
geo = fullfile(repo, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody', 'mjc', 'gait2392_robot', 'Geometry');
lib = fullfile(repo, 'Solid_Models', 'Simscape_Part_Library', 'OpenSim_Bones');
if ~isfolder(lib), mkdir(lib); end
bones = {'femur_r.stl','femur_l.stl','tibia_r.stl','tibia_l.stl','fibula.stl', ...
    'foot.stl','l_foot.stl','bofoot.stl','l_bofoot.stl','talus.stl','l_talus.stl', ...
    'pelvis.stl','sacrum.stl'};
nCopied = 0;
for k = 1:numel(bones)
    src = fullfile(geo, bones{k});
    if isfile(src)
        copyfile(src, fullfile(lib, bones{k}), 'f');
        nCopied = nCopied + 1;
    end
end
fprintf('part library: %d bone STLs in %s\n', nCopied, lib);

% ---- demo model ---------------------------------------------------------------
mdl = 'sns_bone_leg_demo';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
B = @(src, n, pos) add_block(src, [mdl '/' n], 'Position', pos);

B('sm_lib/Frames and Transforms/World Frame', 'World', [30 200 60 230]);
B('nesl_utility/Solver Configuration', 'Solver', [30 300 90 350]);
B('sm_lib/Utilities/Mechanism Configuration', 'Mech', [30 380 120 430]);
set_param([mdl '/Mech'], 'GravityVector', '[0 -9.80665 0]');   % OpenSim y-up

B('sm_lib/Joints/Spherical Joint', 'hip', [160 200 220 260]);
B('sm_lib/Body Elements/File Solid', 'femur_solid', [420 180 520 240]);
B('sm_lib/Joints/Revolute Joint', 'knee', [560 200 620 260]);
B('sm_lib/Body Elements/File Solid', 'tibia_solid', [680 180 780 240]);
B('sm_lib/Joints/Revolute Joint', 'ankle', [820 200 880 260]);
B('sm_lib/Body Elements/File Solid', 'foot_solid', [940 180 1040 240]);

set_bone([mdl '/femur_solid'], fullfile(lib, 'femur_r.stl'));
set_bone([mdl '/tibia_solid'], fullfile(lib, 'tibia_r.stl'));
set_bone([mdl '/foot_solid'], fullfile(lib, 'foot.stl'));

phW = get_param([mdl '/World'], 'PortHandles');
phH = get_param([mdl '/hip'], 'PortHandles'); hP = [phH.RConn phH.LConn];
phF = get_param([mdl '/femur_solid'], 'PortHandles'); fP = [phF.RConn phF.LConn];
phK = get_param([mdl '/knee'], 'PortHandles'); kP = [phK.RConn phK.LConn];
phT = get_param([mdl '/tibia_solid'], 'PortHandles'); tP = [phT.RConn phT.LConn];
phA = get_param([mdl '/ankle'], 'PortHandles'); aP = [phA.RConn phA.LConn];
phS = get_param([mdl '/foot_solid'], 'PortHandles'); sP = [phS.RConn phS.LConn];
phSv = get_param([mdl '/Solver'], 'PortHandles'); svP = [phSv.RConn phSv.LConn];
phM = get_param([mdl '/Mech'], 'PortHandles'); mP = [phM.RConn phM.LConn];

W = phW.RConn(1);
add_line(mdl, svP(1), W);      % branch solver onto the world net
add_line(mdl, mP(1), W);       % branch mechanism config
add_line(mdl, W, hP(1));       % hip base
add_line(mdl, hP(2), fP(1));   % femur on hip follower
add_line(mdl, fP(1), kP(1));   % knee base branches the femur frame
add_line(mdl, kP(2), tP(1));
add_line(mdl, tP(1), aP(1));
add_line(mdl, aP(2), sP(1));

try
    set_param(mdl, 'StopTime', '0.05');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('BONE DEMO: UPDATE OK\n');
    save_system(mdl, fullfile(fileparts(mfilename('fullpath')), '..', 'sns_bone_leg_demo.slx'));
    sim(mdl);
    fprintf('BONE DEMO: SIM OK\n');
catch ME
    fprintf('BONE DEMO FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 6)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 170)));
    end
    save_system(mdl, fullfile(fileparts(mfilename('fullpath')), '..', 'sns_bone_leg_demo.slx'));
end
close_system(mdl, 0);
fprintf('=== sns_bone_pipeline DONE ===\n');
end

function set_bone(blk, stlPath)
% File Solid: geometry + Onyx density (schema-adaptive)
set_param(blk, 'ExtGeomFileName', stlPath);
try, set_param(blk, 'ExtGeomFileUnits', 'm'); catch, end
try, set_param(blk, 'UnitType', 'Custom'); catch
    try, set_param(blk, 'UnitType', 'FromFile'); catch, end
end
try
    set_param(blk, 'InertiaType', 'CalculateFromGeometry');
catch
    try, set_param(blk, 'InertiaType', 'Density'); catch, end
end
try, set_param(blk, 'BasedOnType', 'Density'); catch, end
try, set_param(blk, 'Density', '1200'); catch ME
    fprintf('density set failed on %s: %s\n', blk, ME.message);
    dp = fieldnames(get_param(blk, 'DialogParameters'));
    fprintf('File Solid params: %s\n', strjoin(dp', ' | '));
end
end
