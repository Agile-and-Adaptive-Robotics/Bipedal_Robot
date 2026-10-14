%% Evaluate the run-4 winner under the run-5 model (p4/p5 dead,
%% circumcircle wrap, condyle constraint) and print the margin/strain
%% profile to locate the shortfall.
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

[xFull, ctxCam] = camToExtDesign20mm(xCamR4, ctx, xFixed);
pred = predictKneeExt20mm(xFull, ctxCam);
[cAll, ~] = nonlconExtCam20mm(xCamR4, ctx, xFixed, 200);

humanTgt = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, 'pchip', 'extrap');
marginReq = (abs(pred.TorqueZ(:)) - 1.05*humanTgt(:))./ctx.torqueScale;
marginPlain = (abs(pred.TorqueZ(:)) - humanTgt(:))./humanTgt(:);

fprintf('run-4 design under the new model:\n');
refPts0 = zeros(3,2);   % regularizer reference unused for diagnostics
fprintf('  objective      = %.6g\n', objective_KneeExtCam20mm(xCamR4, ctx, xFixed, refPts0));
fprintf('  max constraint = %+.6f (cRestLength %+.6f, cContraction %+.6f, cVtug %+.6f)\n', ...
    max(cAll), cAll(5), cAll(6), cAll(end));
fprintf('  min margin vs required = %+.2f%%; vs plain = %+.2f%%\n', ...
    100*min(marginReq), 100*min(marginPlain));
fprintf('  min strain_p %.4f, max strain_f %.4f (KMAX %.3f)\n', ...
    min(pred.strain_p), max(pred.strain_f), ctx.KMAX);

[~, iWorst] = min(marginReq);
fprintf('  worst angle %.1f deg: |Tq| %.2f Nm vs required %.2f Nm; strain_f %.3f\n', ...
    ctx.phiD(iWorst), abs(pred.TorqueZ(iWorst)), 1.05*humanTgt(iWorst), pred.strain_f(iWorst));

fprintf('\n  angle    |Tq620|   reqTgt   margin%%   strain_f  arm_mm\n');
for i = 1:5:numel(ctx.phiD)
    fprintf('%7.1f %9.2f %8.2f %8.1f %9.3f %7.1f\n', ctx.phiD(i), ...
        abs(pred.TorqueZ(i)), 1.05*humanTgt(i), 100*marginReq(i), ...
        pred.strain_f(i), 1000*pred.momentArm(i));
end

disp('DONE');
