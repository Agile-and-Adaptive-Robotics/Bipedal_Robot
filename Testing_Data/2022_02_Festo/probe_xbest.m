%probe_xbest.m — verify the load-and-display workflow end to end for BOTH knees:
%rebuild each run's ctx, re-evaluate the saved xBest, confirm the numbers match
%the run's own validation prints. This is the exact state Ben's interactive
%load-and-run-from-the-display-section workflow needs.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization');

%% Extensor: rehydrate from the Vas_Pam capture and verify
V = load('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_20260910_0528.mat', 'xBest', 'XiUsed', 'fBest', 'cBest');
ctxE = buildKneeExtContext20mm();
ctxE.Xi0 = V.XiUsed(1); ctxE.Xi1 = V.XiUsed(2); ctxE.Xi2 = V.XiUsed(3); ctxE.Xi3 = V.XiUsed(4);
predE = predictKneeExt20mm(V.xBest, ctxE);
cE = nonlconExt20mm(V.xBest, ctxE);
JE = objective_KneeExt20mm(V.xBest, ctxE);
fprintf('[Extensor] run fBest %.6g | fresh-ctx re-eval %.6g | max nonlcon %.6g vs %.6g\n', ...
    V.fBest, JE, max(cE), max(V.cBest));
torqueMarginE = (predE.TorqueZ(:) - ctxE.humanTorque(:)) ./ ctxE.humanTorque(:);
fprintf('[Extensor] min torque margin (fresh ctx): %+.3f%% at %.1f deg\n', ...
    100*min(torqueMarginE), ctxE.phiD(abs(torqueMarginE) == min(abs(torqueMarginE))));

%% Flexor: rehydrate from the Bifemsh capture and verify
B = load('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Bifemsh_20mm_Result_20260910_1234.mat', 'xBest', 'fBest', 'cCollision', 'Xi3');
ctxF = buildKneeFlexorContext20mm();
predF = predictKneeFlexor20mm(B.xBest, ctxF);
geo = ctxF.geo; idxP2 = 4:6;
cF = nonlconExclusion(B.xBest, geo, ctxF, idxP2);
JF = objective_KneeFlexor20mm(B.xBest, ctxF);
fprintf('[Flexor] saved fBest %.6g | fresh-ctx re-eval %.6g | max nonlcon saved %.6g fresh %.6g\n', ...
    B.fBest, JF, max(B.cCollision), max(cF));
torqueMarginF = (predF.TorqueZ(:) - ctxF.humanTorque(:)) ./ ctxF.humanTorque(:);
fprintf('[Flexor] min torque margin (fresh ctx): %+.3f%% at %.1f deg\n', ...
    100*min(torqueMarginF), ctxF.phiD(abs(torqueMarginF) == min(abs(torqueMarginF))));
fprintf('PROBE DONE\n');
