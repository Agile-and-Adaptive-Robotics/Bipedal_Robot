%% Verification of the cam design-of-record mat (30-torque drift fix).
%% Read-only: loads Vas_Pam_20mm_Result_cam_20261014_0416.mat and prints
%% its own recorded quantities; reproduces requiredPressureScanCam
%% (copied verbatim from Opt_run_ExtCam.m:376-422) on the saved design.
root = 'D:/Github/Bipedal_Robot';
moDir = fullfile(root,'Code','Matlab','Mesh_Optimization');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(moDir);
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

S = load(fullfile(moDir,'Results','Vas_Pam_20mm_Result_cam_20261014_0416.mat'));
fprintf('mat has xBestCam/xFullBest/ctxCamBest/predBest/margins: %d %d %d %d %d %d\n', ...
    isfield(S,'xBestCam'), isfield(S,'xFullBest'), isfield(S,'ctxCamBest'), ...
    isfield(S,'predBest'), isfield(S,'marginBest'), isfield(S,'marginReqBest'));

mReq = 100*min(S.marginReqBest);
mPlain = 100*min(S.marginBest);
wrapR = 1000*S.predBest.routeInfo.wrap.radius(2,1);
arm = S.predBest.momentArm;
fprintf('FINAL (mat) min margin vs required: %+.2f%%\n', mReq);
fprintf('FINAL (mat) min margin vs plain   : %+.2f%%\n', mPlain);
fprintf('FINAL (mat) wrap bend radius      : %.1f mm\n', wrapR);
fprintf('FINAL (mat) arm min/mean/max      : %.1f / %.1f / %.1f mm\n', ...
    1000*min(arm), 1000*mean(arm), 1000*max(arm));
fprintf('FINAL (mat) p8 (mat coords)       : [%.6f, %.6f], tendon %.5f m\n', ...
    S.xBestCam(5), S.xBestCam(6), S.xBestCam(7));

%% requiredPressureScanCam verbatim reproduction (Opt_run_ExtCam.m:376-422)
xFull = S.xFullBest; ctxEval = S.ctxCamBest;
Pgrid = (0:20:620)';
humanTgt = interp1(ctxEval.humanAngleD, ctxEval.humanTorqueAbs, ...
    ctxEval.phiD, 'pchip', 'extrap');
hasTgt = isfinite(humanTgt) & humanTgt > 0;
restV = xFull(7); kmaxV = restV*(1 - ctxEval.KMAX);
[LocationV, bendV] = buildDistalRingLocation20mm( ...
    xFull(1:3), xFull(4:6), xFull(8), ctxEval);
Tq = nan(numel(Pgrid), numel(ctxEval.phiD));
for k = 1:numel(Pgrid)
    b = MonoPamDataExplicit_balanceX3(ctxEval.Name, LocationV, ...
        ctxEval.CrossPoint, ctxEval.Dia, ctxEval.T_Pam, restV, kmaxV, ...
        xFull(8), ctxEval.fitting, Pgrid(k), ctxEval.Xi0, ctxEval.Xi1, ...
        ctxEval.Xi2, ctxEval.Xi3, ctxEval.wraps, ctxEval.phiD, ...
        ctxEval.BPAcount, bendV);
    Tq(k,:) = abs(b.Torque_p(:,3)).';
end
Preq = nan(numel(ctxEval.phiD),1);
for i = find(hasTgt).'
    col = Tq(:,i);
    fin = isfinite(col);
    if ~any(fin) || max(col(fin)) < humanTgt(i)
        continue
    end
    k = find(col >= humanTgt(i), 1);
    if k <= 1
        Preq(i) = 0;
    else
        pq = Pgrid([k-1 k]);
        [xs, isrt] = sort(col([k-1 k]));
        Preq(i) = interp1(xs, pq(isrt), humanTgt(i));
    end
end
pmin = min(Preq(find(hasTgt)),[],'omitnan');
pmax = max(Preq(find(hasTgt)),[],'omitnan');
fprintf('FINAL (mat) required pressure     : %.0f - %.0f kPa\n', pmin, pmax);
exit
