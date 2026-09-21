%% run_osim_import_and_rig_inventory_20260920.m
% 1) Run the OpenSim gait2392 -> URDF -> smimport route (license-blocked on
%    easteregg2; laptop license carries SimMechanics).
% 2) Inventory the previously imported 09_BA_003 knee-rig model: which joint
%    blocks smimport made from the Multibody-Link XML.

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
addpath(sns);
fprintf('=== MATLAB %s ===\n', version);

%% --- 1. OpenSim -> Simscape Multibody import -------------------------------
try
    run(fullfile(sns, 'sns_osim_import.m'));
catch ME
    fprintf('OSIM IMPORT FAILED: %s\n', ME.message);
end

%% --- 2. knee-rig import inventory -------------------------------------------
rig = fullfile(sns, 'mdl_knee_rig_import_tmp_imported.slx');
try
    load_system(rig);
    mdl = 'mdl_knee_rig_import_tmp';
    blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    fprintf('\nRIG INVENTORY: %d blocks total in %s\n', numel(blks), mdl);

    % joint blocks: smimport names them "<joint> Joint" reference blocks from sm_lib
    joints = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', ...
        'ReferenceBlock', 'sm_lib/Joints/.[jJ]oint.*');
    if isempty(joints)
        % fallback: match by name suffix
        joints = blks(contains(lower({blks.Name}), 'joint'));
    end
    fprintf('joint blocks: %d\n', numel(joints));
    for k = 1:numel(joints)
        rb = '';
        try, rb = get_param(joints{k}, 'ReferenceBlock'); end
        fprintf('  %-40s  <- %s\n', joints{k}, rb);
    end

    % body/solid blocks
    sol = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', ...
        'Name', '.*[Ss]olid.*');
    fprintf('solid-ish blocks: %d\n', numel(sol));
    % 6-DOF / weld / fixed joints
    for pat = {'Weld Joint', '6-DOF Joint', 'Bushing Joint'}
        w = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Name', ['.*' pat{1} '.*']);
        if ~isempty(w), fprintf('%s count: %d\n', pat{1}, numel(w)); end
    end
    close_system(mdl, 0);
    fprintf('RIG INVENTORY DONE\n');
catch ME
    fprintf('RIG INVENTORY FAILED: %s\n', ME.message);
end
fprintf('=== ALL DONE ===\n');
