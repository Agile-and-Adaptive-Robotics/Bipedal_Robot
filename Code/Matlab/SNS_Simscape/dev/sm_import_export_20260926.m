function sm_import_export_20260926(xmlName, mdl)
%% sm_import_export_20260926  Import Ben's fresh Simscape Multibody Link XML
%% exports (2026-09-26) into Simscape Multibody and save to SNS_Simscape.
%%
%%   sm_import_export_20260926('09_BA_003.xml', 'mdl_leg_rig_ba003')
%%   sm_import_export_20260926('10_AH_001.xml', 'mdl_humanoid_lower_ah001')
%%
%% Pattern proven 2026-09-20 (dev\imports2_20260920.m): sanitized same-folder
%% copy (so relative geometry refs resolve), smimport, joint inventory via
%% ReferenceBlock, gravity -Z (CAD Z up), save BEFORE the updateDiagram so a
%% compile failure still leaves the imported model on disk.

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
addpath(sns);
fprintf('=== MATLAB %s | %s ===\n', version, mdl);

xml = fullfile(sns, '..', '..', '..', 'Solid_Models', ...
    'Biomimetics_2022-Knee_Test', 'Knee assembly', xmlName);
assert(isfile(xml), 'XML not found: %s', xml);
srcDir = fileparts(xml);
tmp = fullfile(srcDir, [mdl '.xml']);
if bdIsLoaded(mdl), close_system(mdl, 0); end
copyfile(xml, tmp, 'f');
clob = onCleanup(@() delete(tmp));

% --- import (model name comes from the sanitized file name) ----------------
smimport(tmp);
assert(bdIsLoaded(mdl), 'smimport did not create %s', mdl);

blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
fprintf('\n%s: %d blocks\n', mdl, numel(blks));
nJ = 0; nSol = 0; nFS = 0;
for k = 1:numel(blks)
    rb = '';
    try, rb = get_param(blks{k}, 'ReferenceBlock'); if ~ischar(rb), rb = ''; end, catch, end
    if ~isempty(strfind(rb, 'sm_lib/Joints/'))
        nJ = nJ + 1;
        fprintf('  JOINT  %-58s <- %s\n', strrep(blks{k}, [mdl '/'], ''), rb);
    elseif ~isempty(strfind(rb, 'sm_lib/Body/Elements/File Solid'))
        nFS = nFS + 1;
    elseif ~isempty(strfind(rb, 'sm_lib/Body/Elements/'))
        nSol = nSol + 1;
    end
end
fprintf('inventory: %d joint blocks, %d File Solids, %d other body elements\n', nJ, nFS, nSol);

% sample-check one File Solid geometry path
fs = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'MaskType', 'File Solid');
if ~isempty(fs)
    fn = '';
    try, fn = get_param(fs{1}, 'FileName'); catch, end
    fprintf('file solid sample: %s\n  exists: %d\n', fn, isfile(fn));
end

% --- gravity: CAD Z is up ---------------------------------------------------
try
    mc = find_system(mdl, 'LookUnderMasks', 'all', 'MaskType', 'Mechanism Configuration');
    if ~isempty(mc)
        set_param(mc{1}, 'GravityVector', '[0 0 -9.80665]');
        fprintf('gravity -Z set\n');
    else
        fprintf('no Mechanism Configuration block found\n');
    end
catch ME
    fprintf('gravity tweak skipped: %s\n', ME.message);
end

% --- save BEFORE compile (crash-hazard ordering) ----------------------------
out = fullfile(sns, [mdl '_imported.slx']);
save_system(mdl, out);   % renames model to <mdl>_imported
mdlS = [mdl '_imported'];

% --- compile check -----------------------------------------------------------
try
    set_param(mdlS, 'StopTime', '0.01');
    set_param(mdlS, 'SimulationCommand', 'update');
    fprintf('%s: UPDATE (compile) OK\n', mdlS);
catch ME
    fprintf('%s: UPDATE FAILED: %s\n', mdlS, ME.message);
    for c = 1:numel(ME.cause)
        m = ME.cause{c}.message;
        fprintf('  CAUSE %d: %s\n', c, m(1:min(end, 200)));
    end
end
save_system(mdlS);
close_system(mdlS, 0);
fprintf('=== %s DONE ===\n', mdl);
end
