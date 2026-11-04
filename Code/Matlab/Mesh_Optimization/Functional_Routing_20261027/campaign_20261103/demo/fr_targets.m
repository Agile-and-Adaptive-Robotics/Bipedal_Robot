function target = fr_targets(C, ~, cs)
%FR_TARGETS Human target torque curve for one group (v2, 2026-11-04).
%
% v2 (Ben directive item 2): the target of record is the SIGNED, PER-POSE
% torque curve built by the targets_v2 pipeline from the STOCK human model
% (targets_v2/<group>_signed_targets.json via fr_loadSignedTargets): per
% grid cell, tau = SUM over constituent muscles of F_m(q) * r_m(q), with
% r = OpenSim computeMomentArm (== -dL/dq, verified) and F_m = Fmax *
% active force-length at the rigid-tendon CE length - NOT a constant-Fmax
% moment-arm curve, not a constant torque, and not a sum of separate
% peaks. The loader validates cell count and swept coordinate against the
% case; documented fallbacks cover demo (gate) cases and dense
% verification grids. There is NO silent fallback to the stage-1
% constant-Fmax construction.
%
% Sign convention (targets_v2 README): positive tau drives the primary
% coordinate in its positive direction per the model; knee_angle_r
% positive = extension, i.e. flexion NEGATIVE (this repo's convention).
% (The second argument is the kinematics handle, kept for call-site
% compatibility with fr_finishCase; the loader consumes the case struct.)

% EXPORT MODE (fr_export_demo_grid only): placeholder zeros so the demo
% case build can run before the demo target JSONs exist. Never used by
% the campaign or the gate: the exporter unsets the flag before any
% scoring run, and the placeholder is loudly marked.
if strcmp(getenv('FR_TARGETS_SKIP'), '1')
    nG0 = numel(cs.passes{1}.qStructs);
    target = struct();
    target.tau = zeros(nG0, 1);
    target.qList = cs.passes{1}.qStructs;
    target.peakAbs = 0;
    target.contractMifN = NaN;
    target.mapId = cs.groupMapId;
    target.peakVsContractPct = NaN;
    target.source = 'EXPORT-PLACEHOLDER (fr_export_demo_grid mode)';
    target.tauIsResampled = false;
    return
end

tgtDir = fullfile(C.camp, 'targets_v2');
try
    tau = fr_loadSignedTargets(cs.name, cs);
    srcNote = sprintf('targets_v2 exact (%s)', cs.name);
catch loadErr
    if strcmp(loadErr.identifier, 'fr_loadSignedTargets:Count')
        % DENSE VERIFICATION grid (override gridN): the 9-cell signed
        % target curve of record is resampled onto the dense sweep by
        % coordinate-value interpolation - the same linear convention
        % fr_objective uses for any length mismatch. Documented as
        % verification-mode resampling, never used for scoring runs.
        T = jsondecode(fileread(fullfile(tgtDir, ...
            sprintf('%s_signed_targets.json', cs.name))));
        req = denseReq(T);
        vals = reshape(cs.passes{1}.th1, 1, []);
        tau = interp1(req, T.tau_vec_Nm(:), vals(:), 'linear', 'extrap');
        srcNote = sprintf(['targets_v2 RESAMPLED onto the DENSE ' ...
            'verification grid (%d cells, %s)'], numel(vals), cs.name);
    elseif ~strcmp(loadErr.identifier, 'fr_loadSignedTargets:Missing')
        rethrow(loadErr);   % dof/grid mismatches stay FATAL
    else
        % DEMO-case fallback (gate demonstrator names like hip_flex_ext /
        % sar_layering): reuse the OWNING contract group's signed curve,
        % resampled onto this case's swept grid BY COORDINATE VALUE.
        % Valid only when the swept dof matches the group's primary dof;
        % anything else errors. Logged loudly - never silent.
        tc = jsondecode(fileread(fullfile(C.targetsDir, ...
            'target_contract.json')));
        grpIdx = find([tc.groups.map_id] == cs.groupMapId, 1);
        assert(~isempty(grpIdx), 'fr_targets:noGroup', ...
            'demo case %s (map id %d) maps to no contract group', ...
            cs.name, cs.groupMapId);
        grpName = tc.groups(grpIdx).group;
        T = jsondecode(fileread(fullfile(tgtDir, ...
            sprintf('%s_signed_targets.json', grpName))));
        coord = cs.passes{1}.sweptCoord;
        assert(strcmp(coord, T.primaryDof), 'fr_targets:dof', ...
            ['demo case %s sweeps %s but group %s targets %s - a ' ...
            'resampled target would be physically meaningless'], ...
            cs.name, coord, grpName, T.primaryDof);
        req = denseReq(T);
        vals = reshape(cs.passes{1}.th1, 1, []);
        tau = interp1(req, T.tau_vec_Nm(:), vals(:), 'linear', 'extrap');
        srcNote = sprintf(['targets_v2 RESAMPLED from group %s onto the ' ...
            '%s demo grid (same dof %s; coordinate-value interpolation)'], ...
            grpName, cs.name, coord);
    end
end
fprintf('fr_targets: %s -> %s\n', cs.name, srcNote);
target = struct();
target.tau = tau;
target.qList = cs.passes{1}.qStructs;
target.peakAbs = max(abs(tau));
% informational: contract mifN (stock group sum basis)
tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
target.contractMifN = NaN;
target.mapId = cs.groupMapId;
for g = 1:numel(tc.groups)
    if tc.groups(g).map_id == cs.groupMapId
        target.contractMifN = tc.groups(g).mifN;
    end
end
target.peakVsContractPct = 100 * (target.peakAbs / ...
    max(target.contractMifN, eps) - 1);
% provenance marker consumed by the reports
target.source = srcNote;
target.tauIsResampled = contains(srcNote, 'RESAMPLED');
end

% ==========================================================================
function req = denseReq(T)
% per-cell requested sweep values, tolerant of jsondecode cell/struct shape
if iscell(T.cells)
    req = cellfun(@(cc) cc.requested_rad, T.cells);
else
    req = arrayfun(@(cc) cc.requested_rad, T.cells);
end
req = req(:);
end
