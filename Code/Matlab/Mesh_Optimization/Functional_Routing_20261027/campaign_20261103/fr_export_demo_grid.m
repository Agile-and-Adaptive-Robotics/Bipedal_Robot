function fr_export_demo_grid()
%FR_EXPORT_DEMO_GRID Dump the demo (gate) case grids for the signed-target
% pipeline. For every fr_cases demo case: pass-1 swept coordinate(s), the
% stock muscle list, and the pass-1 qStructs VERBATIM (every coordinate
% value per cell, so paired-profile cases like hip_flex_extrot are
% reconstructed exactly). Writes grid_demo_export.json next to
% grid_export.json in the targets scratch folder.
%
% Run headless:
%   matlab -batch "addpath('<camp>'); fr_export_demo_grid"

C = fr_setup();
if isfield(C, 'kinematics'); kin = C.kinematics; else; kin = fr_kinematics(C.osimStock); end
% export mode: placeholder targets while the demo JSONs do not exist yet
setenv('FR_TARGETS_SKIP', '1');
cln = onCleanup(@() setenv('FR_TARGETS_SKIP', ''));
cases = fr_cases(C, kin);

out = struct();
out.machinery = 'fr_cases demo pass-1 qStructs verbatim (all coordinate values per cell)';
out.generated = datetime('now');
out.cases = {};
for c = 1:numel(cases)
    cs = cases(c);
    p1 = cs.passes{1};
    cells = cell(1, numel(p1.qStructs));
    for i = 1:numel(p1.qStructs)
        q = p1.qStructs{i};
        fn = fieldnames(q);
        vals = struct();
        for k = 1:numel(fn)
            vals.(fn{k}) = q.(fn{k});
        end
        cells{i} = struct('index', i, 'q', vals, ...
            'requested_rad', p1.th1(i), ...
            'requested_deg', p1.th1(i) * 180 / pi);
    end
    out.cases{c} = struct('name', cs.name, 'mapId', cs.groupMapId, ...
        'primaryDof', p1.sweptCoord, ...
        'stockMuscles', {cs.stockMuscles}, ...   % cell value inside braces (struct-array trap)
        'nCells', numel(cells), 'cells', {cells});
end

tgtScratch = fullfile('D:/GitHub/Bipedal_Robot/tmp', ...
    'functional_routing_20261103');
fid = fopen(fullfile(tgtScratch, 'grid_demo_export.json'), 'w');
fwrite(fid, jsonencode(out));
fclose(fid);
fprintf('demo grid export: %d cases -> %s\n', numel(cases), ...
    fullfile(tgtScratch, 'grid_demo_export.json'));
for c = 1:numel(out.cases)
    fprintf('  %-16s dof %-16s muscles {%s}\n', out.cases{c}.name, ...
        out.cases{c}.primaryDof, strjoin(out.cases{c}.stockMuscles, ','));
end
end
