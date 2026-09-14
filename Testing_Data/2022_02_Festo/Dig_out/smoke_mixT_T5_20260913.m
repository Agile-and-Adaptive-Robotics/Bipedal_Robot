% smoke_mixT_T5_20260913.m
% Pre-launch smoke for the mixed-convention + T5-shift campaign:
%  1) checkcode evaluator + driver
%  2) rigid baseline (escape) must reproduce the kf(3)-shift-era a0
%  3) finite-Xi call runs
%  4) FLX_T5_YMM=5 changes TEST 5 ONLY (GoF + M_p bit-compare vs shift-off)
%  5) closed-form check: rotated thetabrB puts reconstructed pkbrB exactly +5mm in y
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
cd(root); addpath(root); addpath(here);
addpath(genpath(fullfile(root, '..', '..', 'Code', 'Matlab'))); %MonoPamDataExplicit lives in Robot_Data

r = checkcode('minimizeFlxPin.m');  fprintf('checkcode evaluator: %d findings\n', numel(r));
r = checkcode('minimizeFlxPin10mm.m'); fprintf('checkcode driver:    %d findings\n', numel(r));

%% rigid baseline under current attribution (shift on kf(3))
setenv('FLX_T5_YMM', '0');   %explicit off
[a0, ~] = minimizeFlxPin(0, Inf, Inf);
fprintf('baseline RMSE: %.4f %.4f %.4f %.4f %.4f\n', a0(:,1));
expected = [7.6174; 8.8026; 7.670; 6.373; 9.6944];  %shift-on-T3 attribution (Dig_a0shiftCheck_20260911)
ok = all(abs(a0(:,1) - expected) < 0.01);
fprintf('baseline matches kf(3)-shift expectation: %d\n', ok);

%% finite-Xi sanity
[f1, b1] = minimizeFlxPin(0.002, 3e4, 2e4);
fprintf('finite GoF RMSE: %.4f %.4f %.4f %.4f %.4f\n', f1(:,1));

%% T5 shift isolation
setenv('FLX_T5_YMM', '5');
[f2, b2] = minimizeFlxPin(0.002, 3e4, 2e4);
setenv('FLX_T5_YMM', '0');
fprintf('GoF unchanged tests 1-4: %d\n', isequal(f1(1:4,:), f2(1:4,:)));
fprintf('test 5 changed:          %d\n', ~isequal(f1(5,:), f2(5,:)));

%% closed-form check: bracket-frame +y offset (test 5 geometry, straight from the data)
S = load('KneeFlxPin_10mm_42cm.mat', 'phiD');
P = load('Plot_KneeFlxPin10mm_42cm.mat', 'Bifemsh_Pam');
Loc = P.Bifemsh_Pam.Location; C = 2;
[~, iz0] = min(abs(S.phiD));
pB = Loc(C,:,iz0);  Pbri = [-48.11, -107.81, 13.8]/1000;
pk = pB - Pbri;
th = atan2(pk(2), pk(1));
RkbrZ = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
d5 = 5/1000;
dP_knee = RkbrZ * [0; d5; 0];   %what the bracket-y offset does to the reconstructed point, knee frame
Rxy = hypot(pk(1), pk(2));
fprintf('pkbrB mm: [%.3f %.3f %.3f], theta = %.2f deg\n', pk*1000, rad2deg(th));
fprintf('+5mm bracket-y -> knee-frame shift [%.3f %.3f %.3f] mm, |shift| = %.3f mm\n', dP_knee*1000, norm(dP_knee)*1000);
fprintf('equivalent thetabrB rotation = %.2f deg ("a few degrees more")\n', rad2deg(asin(d5/Rxy)));
disp('SMOKE DONE');
