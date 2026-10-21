% Extract 20mm redesign route points for the supplement routing tables.
% Read-only over the mats of record; prints values for transcription.
% SCRATCH (tmp/paper_bridge), 2026-10-20 supplement builder.
resDir = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results';

%% Extensor 20 mm cam design of record
S = load(fullfile(resDir,'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), ...
    'predBest','ctxCamBest','xBestCam','xFullBest','refPts','seedBase');
ri = S.predBest.routeInfo;
fprintf('=== EXTENSOR 20mm cam (Vas_Pam_20mm_Result_cam_20261014_0416.mat) ===\n');
fprintf('seedOriginal (9x3, ICR/t1 frame rows):\n'); disp(ri.seedOriginal);
fprintf('seed20 (9x3):\n'); disp(ri.seed20);
fprintf('addedPoint (9x3):\n'); disp(ri.addedPoint);
fprintf('eliminatedPoint (9x3):\n'); disp(ri.eliminatedPoint);
fprintf('eliminatedAngleD=%s\n', mat2str(ri.eliminatedAngleD,7));
fprintf('addedAngleD=%s\n', mat2str(ri.addedAngleD,7));
fprintf('activationOrder=%s eliminationOrder=%s\n', mat2str(ri.activationOrder'), mat2str(ri.eliminationOrder'));
fprintf('CrossPoint=%d routeRows=%d\n', ri.CrossPoint, ri.routeRows);
fprintf('geo.seedP6XY=%s seedP7XY=%s seedP8XY=%s\n', ...
    mat2str(S.ctxCamBest.geo.seedP6XY,7), ...
    mat2str(S.ctxCamBest.geo.seedP7XY,7), ...
    mat2str(S.ctxCamBest.geo.seedP8XY,7));
fprintf('xBestCam=%s\n', mat2str(S.xBestCam,7));
fprintf('xFullBest=%s\n', mat2str(S.xFullBest,7));
fprintf('p1=%s pEnd=%s rest=%.6f tendon=%.5f\n', mat2str(S.predBest.p1,6), ...
    mat2str(S.predBest.pEnd,6), S.predBest.rest, S.predBest.tendon);
[~,i45] = min(abs(ri.sweepD - (-45)));
[~,i0]  = min(abs(ri.sweepD - 0));
fprintf('raw route @ %.2f deg:\n', ri.sweepD(i45)); disp(ri.raw(:,:,i45));
fprintf('raw route @ %.2f deg:\n', ri.sweepD(i0));  disp(ri.raw(:,:,i0));

%% Flexor 20 mm design of record
F = load(fullfile(resDir,'Bifemsh_20mm_Result.mat'), 'xBest','pOptimized', ...
    'pInitialWrapped','pChanged','p1Bfinal','p2Bfinal','Xi3','XiUsed');
fprintf('=== FLEXOR 20mm (Bifemsh_20mm_Result.mat) ===\n');
fprintf('xBest=%s\n', mat2str(F.xBest,7));
fprintf('pOptimized=\n'); disp(F.pOptimized);
fprintf('pInitialWrapped=\n'); disp(F.pInitialWrapped);
fprintf('pChanged=\n'); disp(F.pChanged);
fprintf('p1Bfinal=%s\n', mat2str(F.p1Bfinal,7));
fprintf('p2Bfinal=%s\n', mat2str(F.p2Bfinal,7));
if isfield(F,'XiUsed'), fprintf('XiUsed=%s\n', mat2str(F.XiUsed,7)); end
fprintf('Xi3=%s\n', mat2str(F.Xi3,7));
