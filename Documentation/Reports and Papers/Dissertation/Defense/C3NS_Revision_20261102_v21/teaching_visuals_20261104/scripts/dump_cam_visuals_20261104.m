% dump_cam_visuals_20261104.m
% Defense teaching-visual source dump: extract the ADOPTED cam-extensor
% design of record (Vas_Pam_20mm_Result_cam_20261014_0416.mat) into JSON
% for the Blender scene builder. READ-ONLY: no data is modified or re-run;
% optimizer outputs are only loaded and serialized.
%
% Run headless (laptop DESKTOP-5Q16KE9):
%   matlab -batch "run('C:/Users/Ben/Documents/GitHub/Bipedal_Robot/tmp/dump_cam_visuals_20261104.m')"

outPath = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/tmp/cam_visuals_dump_20261104.json';
meshDir = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Defense/C3NS_Revision_20261102_v21/teaching_visuals_20261104/data';

root = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab'));

matFile = fullfile(root, 'Code/Matlab/Mesh_Optimization/Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat');
S = load(matFile);
ctx = S.ctxCamBest;
predBest = S.predBest;
ri = predBest.routeInfo;
phiD = ctx.phiD(:);
N = numel(phiD);
fprintf('Loaded %s: %d poses, phiD [%.2f .. %.2f] deg\n', matFile, N, min(phiD), max(phiD));

D = struct();

% --- top-level design quantities -----------------------------------------
D.source_mat = 'Vas_Pam_20mm_Result_cam_20261014_0416.mat';
D.n_poses = N;
D.phiD = phiD(:).';
if ismatrix(ri.wrap.radius) && numel(ri.wrap.radius) <= 4000
    D.wrap_radius_full = 1000 * ri.wrap.radius;
else
    D.wrap_radius_full = sprintf('<%s %s>', class(ri.wrap.radius), mat2str(size(ri.wrap.radius)));
end
D.wrap_fields = fieldnames(ri.wrap).';
D.eliminatedAngleD = ri.eliminatedAngleD(:).';
D.active_first = ri.active(:, 1).';
D.active_last  = ri.active(:, end).';

% --- route: 9 points x 3 x N (native frames) -----------------------------
D.route_raw = ri.raw;
D.route_active = ri.active;

% --- transforms per pose --------------------------------------------------
D.T_Pam = ctx.T_Pam;
D.T_ICR_t1 = ctx.T_ICR_t1;
if isfield(ctx, 'T_t1_ICR'), D.T_t1_ICR = ctx.T_t1_ICR; end

% --- geometry envelopes ---------------------------------------------------
geoFields = fieldnames(ctx.geo);
D.geo = struct();
for k = 1:numel(geoFields)
    v = ctx.geo.(geoFields{k});
    if isnumeric(v) && numel(v) <= 4000
        D.geo.(geoFields{k}) = v;
    elseif ischar(v)
        D.geo.(geoFields{k}) = v;
    else
        D.geo.(geoFields{k}) = sprintf('<%s %s skipped>', class(v), mat2str(size(v)));
    end
end
fprintf('geo fields: %s\n', strjoin(geoFields.', ', '));

% --- prediction quantities of record --------------------------------------
D.momentArm_m = predBest.momentArm(:).';
D.momentArmVector = predBest.momentArmVector;
D.TorqueZ = predBest.TorqueZ(:).';
D.strain = predBest.strain(:).';
if isfield(predBest, 'humanAngleD'), D.humanAngleD = predBest.humanAngleD(:).'; end
if isfield(predBest, 'humanTorque'), D.humanTorque = predBest.humanTorque(:).'; end
if isfield(ctx, 'humanAngleD'), D.ctx_humanAngleD = ctx.humanAngleD(:).'; end
if isfield(ctx, 'humanTorque'), D.ctx_humanTorque = ctx.humanTorque(:).'; end
if isfield(ctx, 'requiredTorqueMargin'), D.requiredTorqueMargin = ctx.requiredTorqueMargin; end
if isfield(ctx, 'torqueScale'), D.torqueScale = ctx.torqueScale; end
if isfield(ctx, 'Xi0'), D.Xi = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3]; end
D.routeInfo_fields = fieldnames(ri).';
% wrap struct scalar/1D members
wf = fieldnames(ri.wrap);
D.wrap_summary = struct();
for k = 1:numel(wf)
    v = ri.wrap.(wf{k});
    if isnumeric(v) && numel(v) <= 2000
        D.wrap_summary.(wf{k}) = v;
    else
        D.wrap_summary.(wf{k}) = sprintf('<%s %s>', class(v), mat2str(size(v)));
    end
end

% --- registration constants (from appendix_c_bone_background.m) ----------
D.bone_registration.femur_native_condyle_mm = [6.52, -408.00];
D.bone_registration.tibia_native_groove_mm = [6.00, -31.255];
D.bone_registration.femur_offset_used = ctx.geo.femurProfileCenter - [0.00652, -0.40800];
D.bone_registration.tibia_offset_used = ctx.geo.tibiaLowerCenter - [0.00600, -0.031255];
D.bone_registration.tibia_crop_local_y = -0.250;
D.bone_points.femur = fullfile(root, 'Code/Matlab/Bone_Mesh_Plots/Open_Sim_Bone_Geometry', 'Femur_Mesh_Points.txt');
D.bone_points.tibia = fullfile(root, 'Code/Matlab/Bone_Mesh_Plots/Open_Sim_Bone_Geometry', 'Tibia_Mesh_Points.txt');

% moment arm cross-check at a few poses: perpendicular distance from ICR to
% last-active-contact->p9 line, computed here (validation only, not saved
% as of record).
check = struct('phiD', {}, 'mA_eval', {}, 'mA_geom', {}, 'icr_femur', {}, 'p_last', {}, 'p9', {});
for ii = round(linspace(1, N, 9))
    P = ri.raw(:,:,ii);
    T = ctx.T_Pam(:,:,ii) * ctx.T_ICR_t1(:,:,ii);
    ICR = (T \ [0;0;0;1]).';
    j = 8; while j >= 6 && ~ri.active(j, ii), j = j - 1; end
    pQ = P(j,:); pICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), pQ); pF = RowVecTrans(ctx.T_Pam(:,:,ii), pICR);
    p9Q = P(9,:); p9I = RowVecTrans(ctx.T_ICR_t1(:,:,ii), p9Q); p9F = RowVecTrans(ctx.T_Pam(:,:,ii), p9I);
    v = p9F - pF; v = v / norm(v);
    w = ICR(1:3) - pF;
    mA = norm(w - dot(w, v) * v);
    check(end+1) = struct('phiD', phiD(ii), 'mA_eval', D.momentArm_m(ii), ...
        'mA_geom', mA, 'icr_femur', ICR(1:3), 'p_last', pF, 'p9', p9F); %#ok<AGROW>
end
D.momentarm_crosscheck = check;
for k = 1:numel(check)
    fprintf('check phiD=%+7.2f  mA_eval=%.4f m  mA_geom=%.4f m (ICR->p%d->p9)\n', ...
        check(k).phiD, check(k).mA_eval, check(k).mA_geom);
end

% bone cloud extents (sanity, for scene scaling)
fem = readmatrix(D.bone_points.femur); fem = fem(all(isfinite(fem),2),:);
tib = readmatrix(D.bone_points.tibia); tib = tib(all(isfinite(tib),2),:);
D.bone_extents.femur = [min(fem), max(fem)];
D.bone_extents.tibia = [min(tib), max(tib)];
fprintf('femur cloud %s  x[%.4f %.4f] y[%.4f %.4f] z[%.4f %.4f]\n', mat2str(size(fem)), ...
    D.bone_extents.femur([1 4]), D.bone_extents.femur([2 5]), D.bone_extents.femur([3 6]));
fprintf('tibia cloud %s  x[%.4f %.4f] y[%.4f %.4f] z[%.4f %.4f]\n', mat2str(size(tib)), ...
    D.bone_extents.tibia([1 4]), D.bone_extents.tibia([2 5]), D.bone_extents.tibia([3 6]));

fid = fopen(outPath, 'w');
fwrite(fid, jsonencode(D));
fclose(fid);
fprintf('WROTE %s (%d bytes)\n', outPath, dir(outPath).bytes);

% also persist the momentarm crosscheck as human-readable csv
T2 = struct2table(check);
writetable(T2, fullfile(fileparts(outPath), 'cam_visuals_mA_crosscheck_20261104.csv'));
fprintf('DONE_DUMP_CAM_VISUALS\n');
