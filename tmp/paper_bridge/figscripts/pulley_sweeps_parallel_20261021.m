%% Corrected pulley sweeps: parallel-BPA force x tackle travel gain
% Ben's doctrine (2026-10-21): insertion FORCE gain = nBPA/G, insertion
% TRAVEL gain = G. The tackle divides force and multiplies travel; the
% parallel/mirrored BPA pair (as in the paper's knee routes) supplies the
% force. Configurations of record:
%   cfg1  nBPA=1 G=1  recorded baseline (regression anchor)
%   cfg2  nBPA=2 G=1  parallel pair, force x2
%   cfg3  nBPA=2 G=2  DOCTRINE: net torque ~= 1 BPA, travel x2 (force gain 0.5)
%   cfg4  nBPA=2 G=4  force gain 0.25 with the pair, travel x4
%   cfg5  nBPA=4 G=4  net torque ~= 1 BPA, travel x4 (force gain 0.25)
% Records of record untouched; outputs are new *_parallel_20261021 files.
clear; clc;

root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
resDir = fullfile(root, 'Code/Matlab/Mesh_Optimization/Results');
reqDir = fullfile(root, 'tmp/paper_bridge/reqtorque_20261020');

N1 = 25; pressList = [300 450 620];
cfgN   = [1 2 2 2 4];
cfgG   = [1 1 2 4 4];
cfgTag = {'1BPA_G1', '2BPA_G1', '2BPA_G2', '2BPA_G4', '4BPA_G4'};
cfgForceGain = [1 2 1 0.5 1];      % nBPA/G
cfgTravelGain = [1 1 2 4 4];       % G
nC = numel(cfgN);

sites = struct('name', 'soleus_r', 'joint', 'ankle', ...
    'reqFile', fullfile(reqDir, 'OpenSim_ReqTorque_soleus_r_ankle_20261020.txt'), ...
    'reqCol', 3, 'reqSrc', 'stock gait2392_simbody computeMomentArm x fmax (20261020)');
sites(2).name = 'tib_ant_r'; sites(2).joint = 'ankle';
sites(2).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_tib_ant_r_ankle_20261020.txt');
sites(2).reqCol = 3;
sites(2).reqSrc = sites(1).reqSrc;
sites(3).name = 'bifemsh_r'; sites(3).joint = 'knee';
sites(3).reqFile = ''; sites(3).reqCol = 0;
sites(3).reqSrc = 'recorded mono mat ctx.humanTorque (Bifemsh_20mm_Result_pulley_20260926_1652)';
sites(4).name = 'med_gas_r'; sites(4).joint = 'knee';
sites(4).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_med_gas_r_knee_20261020.txt');
sites(4).reqCol = 3; sites(4).reqSrc = sites(1).reqSrc;
sites(5).name = 'rect_fem_r'; sites(5).joint = 'hip';
sites(5).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_rect_fem_r_hip_flexion_r_20261020.txt');
sites(5).reqCol = 3; sites(5).reqSrc = sites(1).reqSrc;

Mmono = load(fullfile(resDir, 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'), 'ctx');

summaryRows = {};
for s = 1:numel(sites)
    st = sites(s);
    f = dir(fullfile(resDir, sprintf('BiPulley_opensim_%s_20260926_*.mat', st.name)));
    S = load(fullfile(resDir, f(1).name), 'spec', 'T', 'points', 'finalObj', ...
        'baseXi0', 'baseXi1', 'baseXi2', 'baseWraps', 'basePressure');
    spec = S.spec;
    fprintf('\n===== %s (%s, %s) =====\n', st.name, st.joint, f(1).name);

    % regression gate: nBPA=1, G=1 must reproduce the recorded grid
    objRec = BiPam_pulley(spec.name, S.points, spec.cross, 20, S.T, ...
        S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, 0.025, ...
        S.basePressure, S.baseXi0, S.baseXi1, S.baseXi2, S.baseWraps, ...
        struct('nPulleyBPA', 1, 'tackleLineParts', [1 1], 'routingMode', 'moving_via'));
    dmax = 0;
    for k = 1:2
        dmax = max(dmax, max(abs(objRec.Torque_ins(:, :, :, k) - ...
            S.finalObj.Torque_ins(:, :, :, k)), [], 'omitnan'));
    end
    fprintf('REGRESSION nBPA=1 G=1 vs recorded: max |dTorque_ins| = %.3e N*m\n', dmax);

    % primary-joint slice grid (secondary pinned at 0 deg)
    kinslice = spec.kin;
    kinslice.thetaRangesDeg(2, :) = [0 0];
    [Tm, th1m] = buildBiTransforms(struct('kin', kinslice), N1, 2);
    th1m = th1m(:)' * 180 / pi;
    ax1 = spec.kin.axes(1, :) / norm(spec.kin.axes(1, :));

    if st.reqCol > 0
        R = readmatrix(st.reqFile);
        reqAng = R(:, 2); reqTau = R(:, st.reqCol);
    else
        reqAng = Mmono.ctx.humanAngleD; reqTau = Mmono.ctx.humanTorque;
    end
    tauReqM = interp1(reqAng, reqTau, th1m(:), 'linear', 'extrap');

    % sweeps over configurations x pressures
    tauSlice = nan(N1, nC, numel(pressList));
    infeasS = false(N1, nC, numel(pressList));
    slackS  = false(N1, nC, numel(pressList));
    for g = 1:nC
        for pp = 1:numel(pressList)
            om = BiPam_pulley(spec.name, S.points, spec.cross, 20, Tm, ...
                S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, ...
                0.025, pressList(pp), S.baseXi0, S.baseXi1, S.baseXi2, ...
                S.baseWraps, struct('nPulleyBPA', cfgN(g), ...
                'tackleLineParts', [cfgG(g) 1], 'routingMode', 'moving_via'));
            tv = squeeze(om.Torque_ins(:, :, 1, 1));
            tauSlice(:, g, pp) = sum(tv .* ax1, 2);
            infeasS(:, g, pp) = om.PulleyInfeasible(:, 1);
            slackS(:, g, pp) = om.PulleySlack(:, 1);
        end
    end

    % margins + summary
    marginS = nan(N1, nC, numel(pressList));
    for g = 1:nC
        for pp = 1:numel(pressList)
            tr = tauSlice(:, g, pp);
            m = (abs(tr) - abs(tauReqM)) ./ abs(tauReqM);
            m(abs(tauReqM) < 0.5) = NaN;
            m(isnan(tr)) = -Inf;
            marginS(:, g, pp) = m;
            pk = max(abs(tr), [], 'omitnan');
            mn = min(m, [], 'omitnan');
            rom = romWidth(th1m, m >= 0);
            [romC, ~] = romContig(th1m, m >= 0);
            summaryRows{end + 1, 1} = st.name;              %#ok<AGROW>
            summaryRows{end, 2} = st.joint;
            summaryRows{end, 3} = cfgTag{g};
            summaryRows{end, 4} = sprintf('nBPA=%d G=%d fGain=%.2f', cfgN(g), cfgG(g), cfgForceGain(g));
            summaryRows{end, 5} = sprintf('%g', pressList(pp));
            summaryRows{end, 6} = sprintf('%.3f', pk);
            summaryRows{end, 7} = sprintf('%.2f', 100 * mn);
            summaryRows{end, 8} = sprintf('%.2f', rom);
            summaryRows{end, 9} = sprintf('%.2f', romC);
            summaryRows{end, 10} = sprintf('%d', sum(infeasS(:, g, pp)));
            summaryRows{end, 11} = sprintf('%d', sum(slackS(:, g, pp)));
            fprintf(['  %s p=%g: peak=%7.2f  minMgn=%8.1f%%  RoM=%6.1f deg ' ...
                '(contig %6.1f)  infeas %d/%d slack %d\n'], cfgTag{g}, ...
                pressList(pp), pk, 100 * mn, rom, romC, ...
                sum(infeasS(:, g, pp)), N1, sum(slackS(:, g, pp)));
        end
    end
    % sanity: parallel pair at G=1 should roughly double the taut peak
    p1 = max(abs(tauSlice(:, 1, end)), [], 'omitnan');
    p2 = max(abs(tauSlice(:, 2, end)), [], 'omitnan');
    p3 = max(abs(tauSlice(:, 3, end)), [], 'omitnan');
    fprintf('SANITY 620 kPa: 1BPA=%.2f  2BPA_G1=%.2f (ratio %.2f)  2BPA_G2=%.2f (ratio %.2f)\n', ...
        p1, p2, p2 / p1, p3, p3 / p1);

    % doctrine torque surface tau(theta, p): nBPA=2, G=2
    pressSurf = [150:50:600, 620];
    tauSurfP = nan(N1, numel(pressSurf));
    for pp = 1:numel(pressSurf)
        os = BiPam_pulley(spec.name, S.points, spec.cross, 20, Tm, ...
            S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, ...
            0.025, pressSurf(pp), S.baseXi0, S.baseXi1, S.baseXi2, ...
            S.baseWraps, struct('nPulleyBPA', 2, ...
            'tackleLineParts', [2 1], 'routingMode', 'moving_via'));
        tv = squeeze(os.Torque_ins(:, :, 1, 1));
        tauSurfP(:, pp) = sum(tv .* ax1, 2);
    end

    out = struct();
    out.site = st.name; out.joint = st.joint; out.sourceMat = f(1).name;
    out.cfgN = cfgN; out.cfgG = cfgG; out.cfgTag = cfgTag;
    out.cfgForceGain = cfgForceGain; out.cfgTravelGain = cfgTravelGain;
    out.pressList = pressList;
    out.theta1Slice_deg = th1m;
    out.tauSlice = tauSlice; out.infeasSlice = infeasS; out.slackSlice = slackS;
    out.tauReqSlice = tauReqM; out.reqAngle = reqAng; out.reqTau = reqTau;
    out.marginSlice = marginS;
    out.pressSurf = pressSurf; out.tauSurfP = tauSurfP;   % nBPA=2, G=2
    out.regressionMaxdT = dmax;
    matName = fullfile(resDir, sprintf('PulleyParallel_%s_20261021.mat', st.name));
    save(matName, '-struct', 'out');
    fprintf('saved %s\n', matName);
end

csvName = fullfile(resDir, 'PulleyParallel_summary_20261021.csv');
fid = fopen(csvName, 'w');
fprintf(fid, 'site,joint,config,forceGain_note,pressure_kPa,peakAbsTau_Nm,minMargin_pct,romWidth_deg,romContig_deg,nInfeasible_slice,nSlack_slice\n');
for i = 1:size(summaryRows, 1)
    fprintf(fid, '%s\n', strjoin(string(summaryRows(i, :)), ','));
end
fclose(fid);
fprintf('\nsaved %s\nDONE\n', csvName);

%% =================== local functions =========================================
function [T, theta1, theta2] = buildBiTransforms(spec, N1, N2)
kin = spec.kin;
r1 = rangeOrDefault(kin, 1, [-30, 90]);
r2 = rangeOrDefault(kin, 2, [-110, 20]);
theta1 = linspace(r1(1), r1(2), N1) * pi / 180;
theta2 = linspace(r2(1), r2(2), N2) * pi / 180;

thK = []; xK = []; yK = [];
if isfield(kin, 'kneeSlot') && ~isempty(kin.kneeSlot) && kin.kneeSlot > 0
    if isfield(kin, 'kneeFit') && ~isempty(kin.kneeFit) && ~isempty(kin.kneeFit.thetaY)
        thK = kin.kneeFit.thetaY;
        xK = interp1(kin.kneeFit.thetaX, kin.kneeFit.x, thK, 'linear', 'extrap');
        yK = kin.kneeFit.y;
    elseif isfield(kin, 'transTheta') && ~isempty(kin.transTheta)
        thK = kin.transTheta;
        xK = interp1(thK, kin.transX, thK, 'linear', 'extrap');
        yK = kin.transY;
    end
end

T = zeros(4, 4, N1, N2);
for i1 = 1:N1
    p1 = kin.pivots(1, :);
    if kin.kneeSlot == 1 && ~isempty(thK)
        p1 = [interp1(thK, xK, theta1(i1), 'linear', 'extrap'), ...
            interp1(thK, yK, theta1(i1), 'linear', 'extrap'), p1(3)];
    end
    T(:, :, i1, 1) = hingeTrans(kin.types{1}, kin.axes(1, :), p1, theta1(i1));
end
for i2 = 1:N2
    p2 = kin.pivots(2, :);
    if kin.kneeSlot == 2 && ~isempty(thK)
        p2 = [interp1(thK, xK, theta2(i2), 'linear', 'extrap'), ...
            interp1(thK, yK, theta2(i2), 'linear', 'extrap'), p2(3)];
    end
    T(:, :, i2, 2) = hingeTrans(kin.types{2}, kin.axes(2, :), p2, theta2(i2));
end
end

function T = hingeTrans(jtype, axis, pivot, theta)
switch jtype
    case 'fixed'
        R = eye(3);
    otherwise
        R = axisRot(axis, theta);
end
T = RpToTrans(R, pivot(:));
end

function R = axisRot(axis, theta)
a = axis(:) / norm(axis);
Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end

function r = rangeOrDefault(kin, slot, def)
if isfield(kin, 'thetaRangesDeg') && size(kin.thetaRangesDeg, 1) >= slot
    r = kin.thetaRangesDeg(slot, :);
else
    r = def;
end
end

function w = romWidth(x, ok)
w = 0;
for i = 1:numel(x) - 1
    if ok(i) && ok(i + 1)
        w = w + (x(i + 1) - x(i));
    elseif ok(i) ~= ok(i + 1)
        w = w + 0.5 * (x(i + 1) - x(i));
    end
end
end

function [wmax, okAny] = romContig(x, ok)
wmax = 0; okAny = any(ok);
i = 1;
while i <= numel(x)
    if ok(i)
        j = i;
        while j < numel(x) && ok(j + 1)
            j = j + 1;
        end
        wmax = max(wmax, x(j) - x(i));
        i = j + 1;
    else
        i = i + 1;
    end
end
end
