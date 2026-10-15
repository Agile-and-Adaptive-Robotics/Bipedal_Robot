%% Probe: run-4 cam point with longer tendons under the new model --
%% does the feasible-with-margin region exist, and which constraint
%% blocks the original tendon?
root = 'D:/GitHub/Bipedal_Robot';
moDir = fullfile(root,'Code','Matlab','Mesh_Optimization');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(moDir);
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

ctx = buildKneeExtContext20mm();
S = load(fullfile(moDir,'Results','Vas_Pam_20mm_Result.mat'),'xBest');
xRef8 = reshape(S.xBest,1,[]);
xFixed = xRef8(1:7);

xCamR4 = [0.112477, 0.014015, 0.125165, -0.020398, ...
          0.065104, -0.020402, 0.083119, 0.0281];

tendonList = [0.083119, 0.095, 0.102, 0.110, 0.117];

for iT = 1:numel(tendonList)
    xCamT = xCamR4; xCamT(7) = tendonList(iT);
    [xFull, ctxCam] = camToExtDesign20mm(xCamT, ctx, xFixed);
    pred = predictKneeExt20mm(xFull, ctxCam);
    [cAll, ~] = nonlconExtCam20mm(xCamT, ctx, xFixed, 200);

    humanTgt = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, 'pchip', 'extrap');
    marginReq = (abs(pred.TorqueZ(:)) - 1.05*humanTgt(:))./ctx.torqueScale;

    fprintf('tendon %.3f: maxC %+.4f | cRest %+.4f cContr %+.4f cVtug %+.4f | p3clip[%+.1f %+.1f] p6Cond %+.1f | minMargin %+.1f%%\n', ...
        tendonList(iT), max(cAll), cAll(5), cAll(6), cAll(end), ...
        cAll(17), cAll(18), cAll(19), 100*min(marginReq));
end

% Also probe a modestly smaller cam at long tendon (run-4 geometry
% scaled back toward the stock upper contact) to test whether the
% torque shortfall at big cams is excursion-driven.
xCamS = xCamR4;
xCamS(1:2) = [0.095, 0.030];         % p6 closer to stock height/x
xCamS(3:4) = [0.105, -0.008];        % keep window + band
for iT = [0.095, 0.105]
    xCamT = xCamS; xCamT(7) = iT;
    [xFull, ctxCam] = camToExtDesign20mm(xCamT, ctx, xFixed);
    pred = predictKneeExt20mm(xFull, ctxCam);
    [cAll, ~] = nonlconExtCam20mm(xCamT, ctx, xFixed, 200);
    humanTgt = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, 'pchip', 'extrap');
    marginReq = (abs(pred.TorqueZ(:)) - 1.05*humanTgt(:))./ctx.torqueScale;
    fprintf('smaller cam tendon %.3f: maxC %+.4f | cRest %+.4f cVtug %+.4f p6Cond %+.1f | minMargin %+.1f%%\n', ...
        iT, max(cAll), cAll(5), cAll(end), cAll(19), 100*min(marginReq));
end

disp('PROBE DONE');
