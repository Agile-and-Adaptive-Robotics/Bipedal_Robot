%% Doctrine post-processing: re-solved pulley mats -> quotable stats + figures
% Reads the three doctrine re-solve batches (2BPA_G2, 2BPA_G4, 4BPA_G4,
% 2026-10-21) plus the recorded gain-1 baseline mats, re-evaluates each on
% the 25-point primary-joint slice at {300,450,620} kPa, computes margin and
% RoM against the same requirements as the 20261020 sweeps, writes
% PulleyDoctrine_summary_20261021.csv + per-site figure-data mats, and
% regenerates the two paper figures in the campaign style.
clear; clc; close all;

root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
resDir = fullfile(root, 'Code/Matlab/Mesh_Optimization/Results');
reqDir = fullfile(root, 'tmp/paper_bridge/reqtorque_20261020');
figOut = fullfile(root, 'Documentation/Reports and Papers', ...
    'Knee_Torque_Test/Figures');
qcDir = fullfile(root, 'tmp/paper_bridge/figscripts');

N1 = 25; pressList = [300 450 620];
% ladder of record: (tag, nBPA, G, file glob). Baseline first.
tags = {'1BPA_G1', '2BPA_G2', '2BPA_G4', '4BPA_G4'};
cfgN = [1 2 2 4];
cfgG = [1 2 4 4];
globPat = {'BiPulley_opensim_%s_20260926_*.mat', ...
    'BiPulley_opensim_%s_2BPA_G2_20261021_*.mat', ...
    'BiPulley_opensim_%s_2BPA_G4_20261021_*.mat', ...
    'BiPulley_opensim_%s_4BPA_G4_20261021_*.mat'};
nC = numel(tags);

sites = struct('name', 'soleus_r', 'joint', 'ankle', ...
    'reqFile', fullfile(reqDir, 'OpenSim_ReqTorque_soleus_r_ankle_20261020.txt'), ...
    'reqCol', 3);
sites(2).name = 'tib_ant_r'; sites(2).joint = 'ankle';
sites(2).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_tib_ant_r_ankle_20261020.txt');
sites(2).reqCol = 3;
sites(3).name = 'bifemsh_r'; sites(3).joint = 'knee';
sites(3).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_bifemsh_r_knee_20261020.txt');
sites(3).reqCol = 3;
sites(4).name = 'med_gas_r'; sites(4).joint = 'knee';
sites(4).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_med_gas_r_knee_20261020.txt');
sites(4).reqCol = 3;
sites(5).name = 'rect_fem_r'; sites(5).joint = 'hip';
sites(5).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_rect_fem_r_hip_flexion_r_20261020.txt');
sites(5).reqCol = 3;

Mmono = load(fullfile(resDir, 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'), 'ctx');

% campaign style
c1 = '#FFD700'; c2 = '#FFB14E'; c3 = '#FA8775'; c4 = '#EA5F94'; c5 = '#CD34B5';
c6 = '#9D02D7'; c7 = '#0000FF'; c8 = '#000000';
c = {c1; c2; c3; c4; c5; c6; c7; c8};
set(groot, 'defaultAxesFontSize', 12, 'defaultAxesFontWeight', 'bold', ...
    'defaultAxesFontName', 'Arial', 'defaultAxesLineWidth', 2, ...
    'defaultTextFontName', 'Arial', 'defaultTextFontSize', 12, ...
    'defaultTextFontWeight', 'bold');

summaryRows = {};
allData = struct();
for s = 1:numel(sites)
    st = sites(s);
    allData(s).name = st.name; allData(s).joint = st.joint;
    fprintf('\n===== %s (%s) =====\n', st.name, st.joint);

    if st.reqCol > 0
        R = readmatrix(st.reqFile);
        reqAng = R(:, 2); reqTau = R(:, st.reqCol);
    else
        reqAng = Mmono.ctx.humanAngleD; reqTau = Mmono.ctx.humanTorque;
    end

    for g = 1:nC
        f = dir(fullfile(resDir, sprintf(globPat{g}, st.name)));
        if isempty(f)
            fprintf('  %s: NO MAT -- skipped\n', tags{g});
            continue;
        end
        S = load(fullfile(resDir, f(1).name), 'spec', 'T', 'points', ...
            'finalObj', 'baseXi0', 'baseXi1', 'baseXi2', 'baseWraps', 'basePressure');
        spec = S.spec;

        % regression gate: rebuild at this config on the recorded 5x5 grid
        objRec = BiPam_pulley(spec.name, S.points, spec.cross, 20, S.T, ...
            S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, 0.025, ...
            S.basePressure, S.baseXi0, S.baseXi1, S.baseXi2, S.baseWraps, ...
            struct('nPulleyBPA', cfgN(g), ...
            'tackleLineParts', [cfgG(g) 1], 'routingMode', 'moving_via'));
        dmax = 0;
        for k = 1:2
            dmax = max(dmax, max(abs(objRec.Torque_ins(:, :, :, k) - ...
                S.finalObj.Torque_ins(:, :, :, k)), [], 'omitnan'));
        end
        fprintf('  %s REGRESSION vs own finalObj: max |dT| = %.3e\n', tags{g}, dmax);

        kinslice = spec.kin;
        kinslice.thetaRangesDeg(2, :) = [0 0];
        [Tm, th1m] = buildBiTransforms(struct('kin', kinslice), N1, 2);
        th1m = th1m(:)' * 180 / pi;
        ax1 = spec.kin.axes(1, :) / norm(spec.kin.axes(1, :));
        tauReqM = interp1(reqAng, reqTau, th1m(:), 'linear', 'extrap');

        tauSlice = nan(N1, numel(pressList));
        infeasS = false(N1, numel(pressList));
        slackS = false(N1, numel(pressList));
        for pp = 1:numel(pressList)
            om = BiPam_pulley(spec.name, S.points, spec.cross, 20, Tm, ...
                S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, ...
                0.025, pressList(pp), S.baseXi0, S.baseXi1, S.baseXi2, ...
                S.baseWraps, struct('nPulleyBPA', cfgN(g), ...
                'tackleLineParts', [cfgG(g) 1], 'routingMode', 'moving_via'));
            tv = squeeze(om.Torque_ins(:, :, 1, 1));
            tauSlice(:, pp) = sum(tv .* ax1, 2);
            infeasS(:, pp) = om.PulleyInfeasible(:, 1);
            slackS(:, pp) = om.PulleySlack(:, 1);
        end

        % doctrine torque surface at this config for the figure
        pressSurf = [150:50:600, 620];
        tauSurfP = nan(N1, numel(pressSurf));
        for pp = 1:numel(pressSurf)
            os = BiPam_pulley(spec.name, S.points, spec.cross, 20, Tm, ...
                S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, ...
                0.025, pressSurf(pp), S.baseXi0, S.baseXi1, S.baseXi2, ...
                S.baseWraps, struct('nPulleyBPA', cfgN(g), ...
                'tackleLineParts', [cfgG(g) 1], 'routingMode', 'moving_via'));
            tv = squeeze(os.Torque_ins(:, :, 1, 1));
            tauSurfP(:, pp) = sum(tv .* ax1, 2);
        end

        for pp = 1:numel(pressList)
            tr = tauSlice(:, pp);
            m = (abs(tr) - abs(tauReqM)) ./ abs(tauReqM);
            m(abs(tauReqM) < 0.5) = NaN;
            m(isnan(tr)) = -Inf;
            pk = max(abs(tr), [], 'omitnan');
            mn = min(m, [], 'omitnan');
            rom = romWidth(th1m, m >= 0);
            [romC, ~] = romContig(th1m, m >= 0);
            summaryRows{end + 1, 1} = st.name;          %#ok<AGROW>
            summaryRows{end, 2} = st.joint;
            summaryRows{end, 3} = tags{g};
            summaryRows{end, 4} = sprintf('nBPA=%d G=%d fGain=%.2f tGain=%d', ...
                cfgN(g), cfgG(g), cfgN(g) / cfgG(g), cfgG(g));
            summaryRows{end, 5} = sprintf('%g', pressList(pp));
            summaryRows{end, 6} = sprintf('%.3f', pk);
            summaryRows{end, 7} = sprintf('%.2f', 100 * mn);
            summaryRows{end, 8} = sprintf('%.2f', rom);
            summaryRows{end, 9} = sprintf('%.2f', romC);
            summaryRows{end, 10} = sprintf('%d', sum(infeasS(:, pp)));
            summaryRows{end, 11} = sprintf('%d', sum(slackS(:, pp)));
            fprintf(['  %s p=%g: peak=%7.2f  minMgn=%8.1f%%  RoM=%6.1f deg ' ...
                '(contig %6.1f)  infeas %d/%d slack %d\n'], tags{g}, ...
                pressList(pp), pk, 100 * mn, rom, romC, ...
                sum(infeasS(:, pp)), N1, sum(slackS(:, pp)));
        end

        allData(s).(['T_' tags{g}]) = struct('th1m', th1m, 'tauSlice', tauSlice, ...
            'tauSurfP', tauSurfP, 'pressSurf', pressSurf, 'tauReqM', tauReqM, ...
            'regressionMaxdT', dmax, 'sizes', struct('RestingL', ...
            S.finalObj.RestingL, 'Kmax', S.finalObj.Kmax, ...
            'TendonL', S.finalObj.TendonL));
    end
end

csvName = fullfile(resDir, 'PulleyDoctrine_summary_20261021.csv');
fid = fopen(csvName, 'w');
fprintf(fid, 'site,joint,config,forceGain_note,pressure_kPa,peakAbsTau_Nm,minMargin_pct,romWidth_deg,romContig_deg,nInfeasible_slice,nSlack_slice\n');
for i = 1:size(summaryRows, 1)
    fprintf(fid, '%s\n', strjoin(string(summaryRows(i, :)), ','));
end
fclose(fid);
fprintf('\nsaved %s\n', csvName);
save(fullfile(resDir, 'PulleyDoctrine_figdata_20261021.mat'), ...
    'allData', 'tags', 'cfgN', 'cfgG', 'pressList', '-v7.3');
fprintf('saved PulleyDoctrine_figdata_20261021.mat\n');

%% =================== FIGURES (campaign style) ================================
siteTags = {'(A)', '(B)', '(C)', '(D)', '(E)'};
siteNames = {'Soleus', 'Tibialis anterior', 'BFSH', 'Medial gastrocnemius', 'Rectus femoris'};

% ---- Figure 1: torque surfaces tau(theta, p) at the doctrine config 2BPA_G2
figT = figure('Color', 'w', 'Position', [50 50 1250 900]);
tT = tiledlayout(figT, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:numel(sites)
    ax = nexttile(tT, s);
    d = allData(s);
    if isfield(d, 'T_2BPA_G2')
        Zt = abs(d.('T_2BPA_G2').tauSurfP).';
        surf(ax, d.('T_2BPA_G2').th1m, d.('T_2BPA_G2').pressSurf, Zt, ...
            'EdgeColor', 'none');
        hold(ax, 'on');
        % requirement-magnitude iso-contour where shapes overlap
        rp = interp1(d.('T_2BPA_G2').th1m, abs(d.('T_2BPA_G2').tauReqM), ...
            d.('T_2BPA_G2').th1m, 'linear', 'extrap');
        contour3(ax, d.('T_2BPA_G2').th1m, d.('T_2BPA_G2').pressSurf, ...
            repmat(rp(:), 1, numel(d.('T_2BPA_G2').pressSurf)).', [0 0], ...
            'LineStyle', '--', 'LineWidth', 1.5, 'Color', c8);
        hold(ax, 'off');
        title(ax, ['\bf ' siteNames{s}], 'Interpreter', 'tex');
        xlabel(ax, '\bf \theta, \circ', 'Interpreter', 'tex');
        ylabel(ax, '\bf p, kPa', 'Interpreter', 'tex');
        zlabel(ax, '\bf |\tau|, N\cdotm', 'Interpreter', 'tex');
        view(ax, 135, 45);
    end
end
nexttile(tT, 6); axis off;
annotation(figT, 'textbox', [0.72 0.06 0.24 0.14], 'String', ...
    {'\bf 2 BPAs, tackle G = 2', 'force gain 1, travel x2', ...
    'dashed: |\tau| = requirement'}, 'Interpreter', 'tex', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none');
exportgraphics(figT, fullfile(figOut, 'PulleyTorqueSurfaces_20261020.eps'), 'ContentType', 'vector');
exportgraphics(figT, fullfile(qcDir, 'QC_PulleyTorqueSurfaces_doctrine.png'), 'Resolution', 150);
close(figT);

% ---- Figure 2: margin surfaces over (theta x travel-gain ladder) at 620 kPa
figM = figure('Color', 'w', 'Position', [50 50 1250 900]);
tM = tiledlayout(figM, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ladderTags = {'1BPA_G1', '2BPA_G2', '4BPA_G4'};
ladderLbl = {'1 BPA, G=1', '2 BPAs, G=2', '4 BPAs, G=4'};
for s = 1:numel(sites)
    ax = nexttile(tM, s);
    d = allData(s);
    have = cellfun(@(t) isfield(d, ['T_' t]), ladderTags);
    if ~all(have), continue; end
    M = nan(N1, numel(ladderTags));
    for g = 1:numel(ladderTags)
        dd = d.(['T_' ladderTags{g}]);
        tr = dd.tauSlice(:, pressList == 620);
        m = (abs(tr) - abs(dd.tauReqM)) ./ abs(dd.tauReqM);
        m(abs(dd.tauReqM) < 0.5) = NaN;
        m(isnan(tr)) = -Inf;
        M(:, g) = m;
    end
    surf(ax, d.('T_1BPA_G1').th1m, 1:numel(ladderTags), M.', 'EdgeColor', 'none');
    hold(ax, 'on');
    contour3(ax, d.('T_1BPA_G1').th1m, 1:numel(ladderTags), M.', [0 0], ...
        'LineStyle', '--', 'LineWidth', 2, 'Color', c8);
    hold(ax, 'off');
    set(ax, 'YTick', 1:numel(ladderTags), 'YTickLabel', ladderLbl);
    title(ax, ['\bf ' siteNames{s}], 'Interpreter', 'tex');
    xlabel(ax, '\bf \theta, \circ', 'Interpreter', 'tex');
    ylabel(ax, '\bf Transmission', 'Interpreter', 'tex');
    zlabel(ax, '\bf Margin', 'Interpreter', 'tex');
    view(ax, 135, 45);
end
nexttile(tM, 6); axis off;
exportgraphics(figM, fullfile(figOut, 'PulleyMarginSurfaces_20261020.eps'), 'ContentType', 'vector');
exportgraphics(figM, fullfile(qcDir, 'QC_PulleyMarginSurfaces_doctrine.png'), 'Resolution', 150);
close(figM);
fprintf('FIGURES WRITTEN\nDONE\n');

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
