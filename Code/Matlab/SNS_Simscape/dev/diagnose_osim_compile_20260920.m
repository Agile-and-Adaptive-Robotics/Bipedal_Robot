%% diagnose_osim_compile_20260920.m — full error report for the imported
%% Gait2392 Simscape model compile, + correct-name knee-rig inventory.

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
addpath(sns);
fprintf('=== MATLAB %s ===\n', version);

%% --- 1. full compile error of the imported OpenSim model --------------------
mdl = 'Gait2392_simbody_simscape';
try
    load_system(fullfile(sns, 'osim_import', [mdl '.slx']));
    assert(strcmp(get_param(mdl, 'Name'), mdl), ...
        'model name is %s, expected %s', get_param(mdl, 'Name'), mdl);
    set_param(mdl, 'StopTime', '0.01');
    try
        set_param(mdl, 'SimulationCommand', 'update');
        fprintf('OSIM COMPILE OK\n');
    catch ME
        fprintf('OSIM COMPILE FAILED — full report:\n');
        fprintf('%s\n', ME.getReport());
        for c = 1:numel(ME.cause)
            fprintf('CAUSE %d: %s\n', c, ME.cause{c}.message);
        end
    end
    % inventory regardless: joints + links
    joints = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', ...
        'ReferenceBlock', 'sm_lib/Joints/.*');
    fprintf('osim model joint blocks: %d\n', numel(joints));
    close_system(mdl, 0);
catch ME
    fprintf('OSIM DIAG OUTER FAIL: %s\n', ME.message);
end

%% --- 2. knee-rig inventory (correct model name = file basename) -------------
try
    rigMdl = 'mdl_knee_rig_import_tmp_imported';
    load_system(fullfile(sns, [rigMdl '.slx']));
    blks = find_system(rigMdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    fprintf('\nRIG INVENTORY: %d blocks total in %s\n', numel(blks), rigMdl);
    joints = find_system(rigMdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', ...
        'ReferenceBlock', 'sm_lib/Joints/.*');
    if isempty(joints)
        joints = blks(contains(lower({blks.Name}), 'joint'));
    end
    fprintf('joint blocks: %d\n', numel(joints));
    for k = 1:numel(joints)
        rb = '';
        try, rb = get_param(joints{k}, 'ReferenceBlock'); end
        fprintf('  %-45s  <- %s\n', strrep(joints{k}, [rigMdl '/'], ''), rb);
    end
    close_system(rigMdl, 0);
    fprintf('RIG INVENTORY DONE\n');
catch ME
    fprintf('RIG INVENTORY FAILED: %s\n', ME.message);
end
fprintf('=== ALL DONE ===\n');
