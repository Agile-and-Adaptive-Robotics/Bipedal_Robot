%% osim_patch_v2_20260920.m — probe the File Solid mesh parameter name, patch
%% all 19 Visual blocks to absolute STL paths, compile, sim, save.

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
addpath(sns);
fprintf('=== MATLAB %s ===\n', version);

mdl = 'Gait2392_simbody_simscape';
load_system(fullfile(sns, 'osim_import', [mdl '.slx']));
geoDir = fullfile(sns, 'osim_import', 'Geometry');

vis = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', ...
    'MaskType', 'File Solid');
fprintf('found %d File Solids\n', numel(vis));

% probe one block's dialog parameters
b = vis{1};
dp = get_param(b, 'DialogParameters');
fn = fieldnames(dp);
fprintf('DialogParameters of %s:\n', b);
for k = 1:numel(fn)
    v = '';
    try, v = get_param(b, fn{k}); if ~ischar(v), v = mat2str(v); end, catch, end
    fprintf('  %-22s = %s\n', fn{k}, v(1:min(end, 90)));
end

% mesh-path parameter (probed 2026-09-20): 'ExtGeomFileName' = 'Geometry/x.stl'
meshParam = 'ExtGeomFileName';

nPatched = 0;
for k = 1:numel(vis)
    fnv = get_param(vis{k}, meshParam);
    [~, nm, xt] = fileparts(fnv);
    absF = fullfile(geoDir, [nm xt]);
    if isfile(absF)
        set_param(vis{k}, meshParam, absF);
        nPatched = nPatched + 1;
    else
        fprintf('  UNRESOLVED %s (param was "%s")\n', vis{k}, fnv);
    end
end
fprintf('patched %d/%d File Solids\n', nPatched, numel(vis));

set_param(mdl, 'StopTime', '0.01');
set_param(mdl, 'SimulationCommand', 'update');
fprintf('OSIM COMPILE OK\n');

nJ = 0; jt = {};
blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
for k = 1:numel(blks)
    rb = '';
    try, rb = get_param(blks{k}, 'ReferenceBlock'); if ~ischar(rb), rb = ''; end, catch, end
    if ~isempty(strfind(rb, 'sm_lib/Joints/'))
        nJ = nJ + 1; jt{end+1} = rb; %#ok<SAGROW>
    end
end
[uJ, ~, iu] = unique(jt);
for u = 1:numel(uJ), fprintf('  %-42s x%d\n', uJ{u}, sum(iu == u)); end
fprintf('osim model: %d blocks, %d joint refs\n', numel(blks), nJ);

sim(mdl);
fprintf('OSIM 0.01 s SIM OK\n');
save_system(mdl);
close_system(mdl, 0);
fprintf('OSIM MODEL SAVED\n=== DONE ===\n');
