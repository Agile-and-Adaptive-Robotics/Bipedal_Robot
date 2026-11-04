function data = fr_export_grid(outPath)
%FR_EXPORT_GRID Export the campaign per-group pass-1 pose grids to JSON.
%
% Directive item 2 (signed per-pose human torque targets), step 1. The
% consumer (build_signed_targets.py) computes tau_target per grid CELL, so
% the exported cells must match fr_cases EXACTLY - same order, same
% clamping. This function therefore does NOT rebuild the grids itself: it
% CALLS fr_cases27(C, kin) (the campaign's own case builder, which the
% driver and the gate use) and dumps each case's pass-1 qStructs verbatim.
%
% Per group the JSON records:
%   group / map_id / region / stockMuscles / diaStr / nBPA seed
%   primaryDof      = passes{1}.sweptCoord (the dof the primary target
%                     torque acts about; only this coordinate moves)
%   pass1Name       = passes{1}.name
%   secondaryDof    = passes{2}.sweptCoord when a coverage pass exists
%                     (informational; only pass 1 is scored)
%   cells(1..9): index, coord, requested_rad/deg (the fr_cases qStructs
%     value, bit-identical), realized_rad/deg (the value the OpenSim
%     evaluation uses after the geometry_20261103 clamping protocol:
%     clamp(requested) into the model's stored Coordinate range), clamped
%     flag, and the model range in radians.
%   coordDefaults: the model-default value of every campaign coordinate
%     (all cells hold every non-swept coordinate at these defaults,
%     exactly as fr_fillQ27 builds the qStructs and as the
%     geometry_20261103 extraction set its poses).
%
% Clamping is applied ONLY into the model's stored range (the
% geometry-README protocol; the model's own range endpoint is kept for the
% knee -120 deg request). The campaign ROM values come from
% actuator_map_27 rom_table_deg, except knee passes use C.kneeRangeDeg
% (fr_cases27 convention, no hyperextension). The loader
% fr_loadSignedTargets validates consumers against the REQUESTED values.
%
% Usage (matlab -batch):
%   fr_export_grid('D:/GitHub/Bipedal_Robot/tmp/functional_routing_20261103/grid_export.json')
% Default output path when omitted is that same tmp file.

if nargin < 1
    outPath = fullfile('D:/GitHub/Bipedal_Robot', 'tmp', ...
        'functional_routing_20261103', 'grid_export.json');
end

C = fr_setup();
kin = fr_kinematics(C.osimStock);
CASES = fr_cases27(C, kin);   % the SAME grid machinery the driver uses

% group index by contract order (fr_cases27 builds in contract order, so
% CASES(i) is contract group i; assert the mapping defensively)
tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
nG = numel(CASES);
assert(nG == numel(tc.groups), ...
    'fr_export_grid:GridCount', ...
    'fr_cases27 returned %d cases but the contract lists %d groups', ...
    nG, numel(tc.groups));

defaults = struct();
for ic = 1:numel(kin.coordNames)
    nm = kin.coordNames{ic};
    defaults.(nm) = kin.qDefault.(nm);
end

groups = cell(1, nG);
for i = 1:nG
    cs = CASES(i);
    assert(strcmp(cs.group, tc.groups(i).group), ...
        'fr_export_grid:Order', ...
        'case %d is %s but contract group %d is %s', i, cs.group, ...
        i, tc.groups(i).group);
    pass1 = cs.passes{1};
    nCells = numel(pass1.th1);
    coord = pass1.sweptCoord;
    hasRange = isKey(kin.coordRange, coord);
    if hasRange
        rg = kin.coordRange(coord);
    else
        rg = [NaN NaN];
    end
    cells = cell(1, nCells);
    anyClamped = false;
    for k = 1:nCells
        req = pass1.th1(k);
        reqDeg = req * 180 / pi;
        [real, clamped] = fr_clampRad(req, rg);
        anyClamped = anyClamped || clamped;
        cells{k} = struct( ...
            'index', k, ...
            'coord', coord, ...
            'requested_rad', req, ...
            'requested_deg', reqDeg, ...
            'realized_rad', real, ...
            'realized_deg', real * 180 / pi, ...
            'clamped', clamped, ...
            'modelRange_rad', rg);
    end
    secondaryDof = '';
    if numel(cs.passes) >= 2
        secondaryDof = cs.passes{2}.sweptCoord;
    end
    groups{i} = struct( ...
        'group', cs.group, ...
        'map_id', cs.groupMapId, ...
        'region', cs.region, ...
        'stockMuscles', {cs.stockMuscles}, ...
        'diaStr', cs.diaStr, ...
        'nBPA_seed', cs.seed.nBPA, ...
        'primaryDof', coord, ...
        'pass1Name', pass1.name, ...
        'secondaryDof', secondaryDof, ...
        'nCells', nCells, ...
        'cells', {cells}, ...
        'coordDefaults_rad', defaults);
    fprintf(['fr_export_grid: %-10s dof=%-16s %d cells ' ...
        '[%.3f .. %.3f] deg clamped=%d\n'], cs.group, coord, nCells, ...
        cells{1}.requested_deg, cells{nCells}.requested_deg, ...
        anyClamped);
end

data = struct( ...
    'contract', 'functional_routing_targets_27_corrected', ...
    'gridMachinery', ['fr_cases27(fr_setup(), fr_kinematics(' ...
    'gait2392_simbody.osim)); pass-1 qStructs verbatim'], ...
    'osimModelStock', C.osimStock, ...
    'targetContract', fullfile(C.targetsDir, 'target_contract.json'), ...
    'actuatorMap', fullfile(C.repo, 'Code', 'Matlab', ...
    'Mesh_Optimization', 'Routing27_v2', 'actuator_map_27.json'), ...
    'kneeRangeOverrideDeg', C.kneeRangeDeg, ...
    'clampingProtocol', ['clamp requested radians into the model''s ' ...
    'stored Coordinate range (geometry_20261103 README protocol; model ' ...
    'clamped flags are all false, so the clamp is applied by the ' ...
    'evaluator, not by OpenSim); only the knee -120 deg request reaches ' ...
    'the stored range endpoint'], ...
    'nGroups', nG, ...
    'groups', {groups});

[fid, msg] = fopen(outPath, 'w');
if fid < 0
    error('fr_export_grid:Open', 'cannot open %s: %s', outPath, msg);
end
fwrite(fid, jsonencode(data, 'PrettyPrint', true));
fclose(fid);
fprintf('fr_export_grid: wrote %s (%d groups)\n', outPath, nG);
end

% ==========================================================================
function [v, clamped] = fr_clampRad(req, rg)
% clamp req into [rg(1) rg(2)] when both endpoints are finite (the
% geometry-README protocol); NaN endpoints (no stored range) pass through.
if any(isnan(rg))
    v = req;
    clamped = false;
    return
end
v = min(max(req, rg(1)), rg(2));
clamped = v ~= req;
end
