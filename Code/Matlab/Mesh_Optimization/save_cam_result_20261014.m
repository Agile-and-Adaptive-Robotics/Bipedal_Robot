%% Mint the dated cam result mat for the FINAL design of record
%% (2026-10-14). The Opt_run_ExtCam runs left their winners in console
%% logs only (save block follows Ben's commented liveRun convention),
%% so this script captures the final design -- run-8 p6/p7 plus the
%% verified p8 clip-edge / +2.05 mm tendon adjustment -- under the
%% p7x-p9x wrap-radius rule, with the same report quantities the
%% driver prints, and saves a full-workspace dated mat in Results.
root = 'D:/GitHub/Bipedal_Robot';
moDir = fullfile(root,'Code','Matlab','Mesh_Optimization');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(moDir);
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

vtugFloorKPa = 150;

ctx = buildKneeExtContext20mm();
S = load(fullfile(moDir,'Results','Vas_Pam_20mm_Result.mat'),'xBest','XiUsed','ctx');
xRef8 = reshape(S.xBest,1,[]);
xFixed = xRef8(1:7);
tendonRef = xRef8(8);

[~,~,routeInfoBase] = buildDistalRingLocation20mm( ...
    xFixed(1:3), xFixed(4:6), tendonRef, ctx);
seedBase = routeInfoBase.seed20;
thetaP6BaseD = atan2d(seedBase(6,2), seedBase(6,1));

adjPush = 0.0254;
refPts = seedBase(6:8,1:2);
refPts(1,:) = refPts(1,:) + adjPush*refPts(1,:)/norm(refPts(1,:));
refPts(2,:) = refPts(2,:) + adjPush*refPts(2,:)/norm(refPts(2,:));

% FINAL design of record: run-8 winner geometry + the verified
% clip-edge p8 (10 mm past the pEnd->p7 chord) and +2.05 mm tendon.
xBestCam = [0.116323, 0.015040, ...
            0.131323, -0.018974, ...
            0.127039, -0.047873, ...
            0.092490];

[xFullBest, ctxCamBest] = camToExtDesign20mm(xBestCam, ctx, xFixed);
predBest = predictKneeExt20mm(xFullBest, ctxCamBest);
[cBest, ~] = nonlconExtCam20mm(xBestCam, ctx, xFixed, vtugFloorKPa);
fBestCam = objective_KneeExtCam20mm(xBestCam, ctx, xFixed, refPts);

humanBest = interp1(ctx.humanAngleD, ctx.humanTorque, ctx.phiD, 'pchip', 'extrap');
marginBest = (predBest.TorqueZ(:) - humanBest(:))./humanBest(:);
requiredBest = (1 + ctx.requiredTorqueMargin)*abs(humanBest(:));
marginReqBest = (abs(predBest.TorqueZ(:)) - requiredBest(:))./ctx.torqueScale;
[stepMax, iStep] = max(abs(diff(predBest.momentArm)));

fprintf('FINAL design under the p7x-p9x wrap-radius rule:\n');
fprintf('wrap bend radius: %.1f mm (p7x-p9x = %.1f mm)\n', ...
    1000*predBest.routeInfo.wrap.radius(2,1), ...
    1000*abs(seedBase(7,1)*0 + xBestCam(3) - xFixed(4)));
fprintf('cRestLength = %+.3f mm; max constraint %+.6f\n', ...
    1000*predBest.cRestLength, max(cBest));
fprintf('min margin vs required %+.2f%%; vs plain %+.2f%%\n', ...
    100*min(marginReqBest), 100*min(marginBest));
fprintf('max arm step %.2f mm at %.1f deg\n', 1000*stepMax, ctx.phiD(iStep));
fprintf('p8 eliminated at: %+.2f (NaN = never)\n', ...
    predBest.routeInfo.eliminatedAngleD(8));

% Full-workspace dated save, Opt_run family pattern. NOT a liveRun mint:
% the design is the already-verified post-hoc final (documented above).
XiUsed = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
stamp = char(string(datetime('now'),'yyyyMMdd_HHmm'));
resultFile = fullfile(moDir,'Results', ...
    sprintf('Vas_Pam_20mm_Result_cam_%s.mat', stamp));
clear root moDir vtugFloorKPa S adjPush iStep stepMax humanBest
save(resultFile);
fprintf('Saved %s\n', resultFile);
