%% osim_recompile_and_inventories_20260920.m — with the Geometry junction in
%% place, recompile the imported OpenSim model; joint inventories for BOTH the
%% osim model and the 09_BA_003 knee rig. Rig inventory runs FIRST; the osim
%% compile LAST (the 2026-09-20 crash in physmod_sm_gui_app_tree.dll happened
%% while closing a model whose update had failed — never touch a model after a
%% failed update in the same process).

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
addpath(sns);
fprintf('=== MATLAB %s ===\n', version);

%% --- 1. knee-rig inventory (09_BA_003 Multibody-Link XML import) ------------
try
    rigMdl = 'mdl_knee_rig_import_tmp_imported';
    load_system(fullfile(sns, [rigMdl '.slx']));
    blks = find_system(rigMdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    fprintf('\nRIG INVENTORY: %d blocks total in %s\n', numel(blks), rigMdl);
    nJ = 0; nSol = 0;
    for k = 1:numel(blks)
        rb = getParamSafe(blks{k}, 'ReferenceBlock');
        if ~isempty(strfind(rb, 'sm_lib/Joints/'))
            nJ = nJ + 1;
            fprintf('  JOINT  %-42s <- %s\n', strrep(blks{k}, [rigMdl '/'], ''), rb);
        elseif ~isempty(strfind(rb, 'sm_lib/Body/')) || ~isempty(strfind(rb, 'sm_lib/Utilities/'))
            nSol = nSol + 1;
        end
    end
    fprintf('rig: %d joint blocks, %d body/utility refs\n', nJ, nSol);
    close_system(rigMdl, 0);
catch ME
    fprintf('RIG INVENTORY FAILED: %s\n', ME.message);
end

%% --- 2. osim model recompile (Geometry junction now resolves the 19 STLs) --
try
    mdl = 'Gait2392_simbody_simscape';
    load_system(fullfile(sns, 'osim_import', [mdl '.slx']));
    set_param(mdl, 'StopTime', '0.01');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('OSIM COMPILE OK (meshes resolved)\n');
    nJ = 0; nPrim = 0; nCont = 0;
    blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    for k = 1:numel(blks)
        rb = getParamSafe(blks{k}, 'ReferenceBlock');
        if ~isempty(strfind(rb, 'sm_lib/Joints/'))
            nJ = nJ + 1;
        elseif ~isempty(strfind(rb, 'sm_lib/Body/'))
            nPrim = nPrim + 1;
        end
    end
    fprintf('osim model: %d blocks, %d joint refs, %d body refs\n', numel(blks), nJ, nPrim);
    % list the joint types compactly
    jt = {};
    for k = 1:numel(blks)
        rb = getParamSafe(blks{k}, 'ReferenceBlock');
        if ~isempty(strfind(rb, 'sm_lib/Joints/'))
            jt{end+1} = rb; %#ok<SAGROW>
        end
    end
    [uJ, ~, iu] = unique(jt);
    for u = 1:numel(uJ)
        fprintf('  %-40s x%d\n', uJ{u}, sum(iu == u));
    end
    % quick 0.01 s sim to prove the plant steps
    out = sim(mdl);
    fprintf('OSIM 0.01 s SIM OK\n');
    close_system(mdl, 0);
catch ME
    fprintf('OSIM RECOMPILE FAILED: %s\n', ME.message);
    for c = 1:numel(ME.cause)
        m = ME.cause{c}.message;
        fprintf('CAUSE %d: %s\n', c, m(1:min(end, 160)));
    end
    % do NOT close or touch the model further (crash hazard after failed update)
end
fprintf('=== ALL DONE ===\n');

function v = getParamSafe(b, p)
v = '';
try
    v = get_param(b, p);
    if ~ischar(v), v = ''; end
catch
end
end
