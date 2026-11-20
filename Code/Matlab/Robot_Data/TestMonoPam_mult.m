function TestMonoPam_mult
%TESTMONOPAM_MULT Test two actual BPA routes with saved knee-result files.
%
% Required on the MATLAB path:
%   MonoPam_mult.m
%   MonoPamDataExplicit_balanceX3.m
%   Vas_Pam_20mm_Result.mat
%   Bifemsh_20mm_Result.mat
%   festo4.m, maxBPAforce.m, RowVecTrans.m, and FestoLookup.mat
%
% Run from the folder containing the files, or add that folder to the path:
%
%   clear classes
%   TestMonoPam_mult

%% Load the saved results
vasFile = which('Vas_Pam_20mm_Result.mat');
flxFile = which('Bifemsh_20mm_Result.mat');
assert(~isempty(vasFile), ...
    'Vas_Pam_20mm_Result.mat was not found on the MATLAB path.')
assert(~isempty(flxFile), ...
    'Bifemsh_20mm_Result.mat was not found on the MATLAB path.')

vas = load(vasFile,'ctx','predBest');
flx = load(flxFile,'ctx','predBest');

%% Extensor check: equivalent parallel paths separated by +/-30 mm in z
v = vas.predBest;
c = vas.ctx;
N = numel(c.phiD);

vasX3 = MonoPamDataExplicit_balanceX3( ...
    'Vastus X3 reference',v.Location,c.CrossPoint,c.Dia,c.T_Pam, ...
    v.rest,v.kmax,v.tendon,c.fitting,c.targetPressure, ...
    c.Xi0,c.Xi1,c.Xi2,c.Xi3,c.wraps,c.phiD,2,v.bendMeasure);

Location1 = v.Location;
Location2 = v.Location;
Location1(:,3,:) = Location1(:,3,:)+0.030;
Location2(:,3,:) = Location2(:,3,:)-0.030;

vasMult = MonoPam_mult( ...
    'Vastus separated pair',{Location1;Location2},c.CrossPoint,c.Dia, ...
    c.T_Pam,v.rest,v.kmax,v.tendon,c.fitting,c.targetPressure, ...
    c.Xi0,c.Xi1,c.Xi2,c.Xi3,c.wraps,c.phiD,2,v.bendMeasure);

assert(isequal(size(vasMult.Torque_p),[N,3]))
assert(all(isfinite(vasMult.Fmag)))
assert(max(abs(vasMult.gama{1}-vasMult.gama{2})) < 1e-12)
assert(max(abs(vecnorm(vasMult.F_p{1},2,2)- ...
    vecnorm(vasMult.F_p{2},2,2))) < 1e-10)

vasTorqueScale = max(abs(vasX3.Torque_p(:,3)),[],'omitnan');
vasTorqueError = max(abs(vasMult.Torque_p(:,3)- ...
    vasX3.Torque_p(:,3)),[],'omitnan')/max(vasTorqueScale,eps);
vasForceMismatch = max(vasMult.forceMismatch./max(vasMult.Fmag,eps), ...
    [],'omitnan');

% The +/-z translations preserve route lengths and sagittal moment arms.
% Report, rather than hide, a substantial difference from the established
% X3 parallel-BPA calculation.
if vasTorqueError >= 0.10
    warning('TestMonoPam_mult:VastusDifference', ...
        'Separated extensor torque differs from X3 by %.3f%%.', ...
        100*vasTorqueError)
end

%% Flexor check: use the saved transform and mirror the route about z = 0
v = flx.predBest;
c = flx.ctx;
N = numel(c.phiD);

flxSingle = MonoPamDataExplicit_balanceX3( ...
    'Bifemsh single-route reference',v.Location,c.CrossPoint,c.Dia, ...
    c.T_Pam,v.rest,v.kmax,v.tendon,c.fitting,c.targetPressure, ...
    c.Xi0,c.Xi1,c.Xi2,c.Xi3,c.wraps,c.phiD,1,v.bendMeasure);

Location1 = v.Location;
Location2 = v.Location;
Location2(:,3,:) = -Location2(:,3,:);

flxMult = MonoPam_mult( ...
    'Bifemsh mirrored pair',{Location1;Location2},c.CrossPoint,c.Dia, ...
    c.T_Pam,v.rest,v.kmax,v.tendon,c.fitting,c.targetPressure, ...
    c.Xi0,c.Xi1,c.Xi2,c.Xi3,c.wraps,c.phiD,2,v.bendMeasure);

assert(isequal(size(flxMult.Torque_p),[N,3]))
assert(all(isfinite(flxMult.Fmag)))
assert(max(abs(flxMult.gama{1}-flxMult.gama{2})) < 1e-12)
assert(max(abs(vecnorm(flxMult.F_p{1},2,2)- ...
    vecnorm(flxMult.F_p{2},2,2))) < 1e-10)

highFlexion = c.phiD(:) <= -90;
if ~any(highFlexion)
    highFlexion = true(N,1);
end

singleTorque = abs(flxSingle.Torque_p(:,3));
pairedTorque = abs(flxMult.Torque_p(:,3));
torqueRatio = pairedTorque./max(singleTorque,eps);
relativeMismatch = flxMult.forceMismatch./max(flxMult.Fmag,eps);
singleForce = vecnorm(flxSingle.F_p,2,2);
pairForceRatio = (2*flxMult.Fmag)./max(singleForce,eps);
singleTorquePerForce = singleTorque./max(singleForce,eps);
pairTorquePerForce = pairedTorque./max(2*flxMult.Fmag,eps);
geometryRatio = pairTorquePerForce./max(singleTorquePerForce,eps);

fprintf('\nTestMonoPam_mult passed its structural checks.\n')
fprintf('Vastus samples                                  = %d\n', ...
    numel(vas.ctx.phiD))
fprintf('Maximum Vastus torque difference from X3       = %.4f %%\n', ...
    100*vasTorqueError)
fprintf('Maximum Vastus relative force mismatch         = %.4f %%\n', ...
    100*vasForceMismatch)
fprintf('Bifemsh samples                                 = %d\n',N)
fprintf('Minimum paired/single torque ratio, <= -90 deg = %.6f\n', ...
    min(torqueRatio(highFlexion),[],'omitnan'))
fprintf('Mean paired/single torque ratio, <= -90 deg    = %.6f\n', ...
    mean(torqueRatio(highFlexion),'omitnan'))
fprintf('Minimum paired/single force ratio, <= -90 deg  = %.6f\n', ...
    min(pairForceRatio(highFlexion),[],'omitnan'))
fprintf('Mean paired/single force ratio, <= -90 deg     = %.6f\n', ...
    mean(pairForceRatio(highFlexion),'omitnan'))
fprintf('Minimum torque-per-force geometry ratio        = %.6f\n', ...
    min(geometryRatio(highFlexion),[],'omitnan'))
fprintf('Mean torque-per-force geometry ratio           = %.6f\n', ...
    mean(geometryRatio(highFlexion),'omitnan'))
fprintf('Maximum Bifemsh relative force mismatch        = %.4f %%\n', ...
    100*max(relativeMismatch,[],'omitnan'))
fprintf('Route 1 strain_f range                         = %.6f to %.6f\n', ...
    min(flxMult.strain_f{1}),max(flxMult.strain_f{1}))
fprintf('Route 2 strain_f range                         = %.6f to %.6f\n', ...
    min(flxMult.strain_f{2}),max(flxMult.strain_f{2}))
fprintf('Minimum route 1 sagittal moment arm            = %.6f m\n', ...
    min(vecnorm(flxMult.mA_p{1}(:,1:2),2,2)))
fprintf('Minimum route 2 sagittal moment arm            = %.6f m\n', ...
    min(vecnorm(flxMult.mA_p{2}(:,1:2),2,2)))
end
