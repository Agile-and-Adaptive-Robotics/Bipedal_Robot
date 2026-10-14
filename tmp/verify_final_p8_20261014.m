%% Verify the final design WITHOUT re-optimizing: run-8 p6/p7, p8 as a
%% 10 mm protrusion past the pEnd->p7 chord (clip distal edge), tendon
%% grown by the added path length. Reports margin, no-stretch, arm
%% continuity, and the p8 contact schedule.
root = 'D:/GitHub/Bipedal_Robot';
moDir = fullfile(root,'Code','Matlab','Mesh_Optimization');
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(moDir);
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

ctx = buildKneeExtContext20mm();
S = load(fullfile(moDir,'Results','Vas_Pam_20mm_Result.mat'),'xBest');
xBest = reshape(S.xBest,1,[]);
p1 = xBest(1:3); pEnd = xBest(4:6); rest = xBest(7);

p6R8 = [0.116323, 0.015040];
p7R8 = [0.131323, -0.018974];
tendonR8 = 0.090440;

u = (p7R8 - pEnd(1:2)); u = u/norm(u);
n = [u(2), -u(1)];
foot = pEnd(1:2) + min(max(dot([0.06836,-0.01599]-pEnd(1:2),u),0), ...
    norm(p7R8-pEnd(1:2)))*u;
p8New = foot + 0.010*n;
tendonNew = tendonR8 + 0.00205;   % + the measured added path length

ctxAdj = ctx;
ctxAdj.geo.seedP6XY = p6R8;
ctxAdj.geo.seedP7XY = p7R8;
ctxAdj.geo.seedP8XY = p8New;
ctxAdj.geo.tibiaWrapCircumcircle = true;

pred = predictKneeExt20mm([p1, pEnd, rest, tendonNew], ctxAdj);
[cAll, ~] = nonlconExtCam20mm([p6R8, p7R8, p8New, tendonNew], ctx, ...
    xBest(1:7), 150);

humanTgt = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, 'pchip', 'extrap');
marginReq = (abs(pred.TorqueZ(:)) - 1.05*humanTgt(:))./ctx.torqueScale;
[stepMax, iStep] = max(abs(diff(pred.momentArm)));

fprintf('p8New t1 = [%.4f %.4f] m, tendon = %.6f m\n', p8New, tendonNew);
fprintf('cRestLength = %+.3f mm (<=0 feasible)\n', 1000*pred.cRestLength);
fprintf('max constraint = %+.6f\n', max(cAll));
fprintf('min margin vs required = %+.2f%%\n', 100*min(marginReq));
fprintf('max arm step %.2f mm at %.1f deg\n', 1000*stepMax, ctx.phiD(iStep));
fprintf('p8: active@fullFlexion %d, eliminated at %+.2f deg (999=never active)\n', ...
    pred.routeInfo.active(8,1), ...
    (pred.routeInfo.eliminatedAngleD(8) + 0*pred.routeInfo.eliminatedAngleD(8)));
fprintf('p8 release angle: %+.2f\n', pred.routeInfo.eliminatedAngleD(8));
fprintf('wrap radius row: %.1f - %.1f mm\n', ...
    1000*min(pred.routeInfo.wrap.radius(2,:)), ...
    1000*max(pred.routeInfo.wrap.radius(2,:)));

disp('VERIFY DONE');
