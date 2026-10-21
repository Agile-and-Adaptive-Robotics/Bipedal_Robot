% pulley_sweeps_20261020.m -- EDIT_SPEC_20261020 section 5 deliverables 1-2:
% fixed-gain EVALUATION sweeps over the recorded 2026-09-26 bi-pulley batch.
% No optimizer runs here; records of record are only READ. New outputs only:
%   Code/Matlab/Mesh_Optimization/Results/PulleySurfaces_<site>_20261020.mat
%   Code/Matlab/Mesh_Optimization/Results/PulleySurfaces_summary_20261020.csv
%
% Per site (med_gas_r, soleus_r, tib_ant_r, rect_fem_r, bifemsh_r):
%   * regression gate: rebuild BiPam_pulley at gain [1 1], 620 kPa on the
%     RECORDED 5x5 T grid with the recorded finalObj sizes; must reproduce
%     the recorded finalObj.Torque_ins (max |delta| reported).
%   * dense grid: T rebuilt inline with the Opt_run_BiPulley two-hinge
%     logic (buildBiTransforms copy, lines 585-634): N1 = 25 on the primary
%     joint, N2 = 5 on the secondary; plus a 25x1 margin slice with the
%     secondary joint pinned at 0 deg.
%   * sweeps: gpair {[1 1],[2 1],[3 1]} x pressure {300,450,620} kPa,
%     routingMode 'moving_via', nPulleyBPA 1 -- the tackle sits at
%     crossing 1 (the live joint of every site).
%   * tau_ax = Torque_ins projected on the hinge axis (signed, N*m);
%     margins vs tau_req (OpenSim required torque; bifemsh = recorded mono
%     ctx.humanTorque); RoM width where margin >= 0 per (G, p).
%
% Bowden note (audit, PULLEY_FACTS section 2): the 09-26 batch used a
% MOVING FORCE POINT (routingMode 'moving_via'); the bowden mode is
% implemented but never exercised on record. These sweeps stay moving_via.

root = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
reqDir = fullfile(root, 'tmp', 'paper_bridge', 'reqtorque_20261020');

gainList = {[1 1], [2 1], [3 1]};
gainTags = {'1/1', '2/1', '3/1'};
pressList = [300 450 620];
N1 = 25;     % dense primary-joint angles
N2 = 5;      % secondary-joint angles

% site table: primary joint label, required-torque source, batch mat
sites = struct('name', 'soleus_r', 'joint', 'ankle', 'cross', 1, ...
    'reqFile', fullfile(reqDir, 'OpenSim_ReqTorque_soleus_r_ankle_20261020.txt'), ...
    'reqCol', 3, 'reqSrc', 'stock gait2392_simbody computeMomentArm x fmax (20261020)');
sites(2).name = 'tib_ant_r';  sites(2).joint = 'ankle'; sites(2).cross = 1;
sites(2).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_tib_ant_r_ankle_20261020.txt');
sites(2).reqCol = 3;
sites(2).reqSrc = 'stock gait2392_simbody computeMomentArm x fmax (20261020)';
sites(3).name = 'bifemsh_r';  sites(3).joint = 'knee';  sites(3).cross = 1;
sites(3).reqFile = '';
sites(3).reqCol = 0;
sites(3).reqSrc = 'recorded mono mat ctx.humanTorque (Bifemsh_20mm_Result_pulley_20260926_1652)';
sites(4).name = 'med_gas_r';  sites(4).joint = 'knee';  sites(4).cross = 1;
sites(4).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_med_gas_r_knee_20261020.txt');
sites(4).reqCol = 3;
sites(4).reqSrc = 'stock gait2392_simbody computeMomentArm x fmax (20261020)';
sites(5).name = 'rect_fem_r'; sites(5).joint = 'hip';   sites(5).cross = 1;
sites(5).reqFile = fullfile(reqDir, 'OpenSim_ReqTorque_rect_fem_r_hip_flexion_r_20261020.txt');
sites(5).reqCol = 3;
sites(5).reqSrc = 'stock gait2392_simbody computeMomentArm x fmax (20261020)';
% med_gas ankle-crossing requirement (secondary crossing), for the record:
reqMedGasAnkle = readmatrix(fullfile(reqDir, ...
    'OpenSim_ReqTorque_med_gas_r_ankle_20261020.txt'));

% mono bifemsh mat = recorded requirement for bifemsh
Mmono = load(fullfile(resDir, 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'), 'ctx');

summaryRows = {};
for s = 1:numel(sites)
    st = sites(s);
    f = dir(fullfile(resDir, sprintf('BiPulley_opensim_%s_*.mat', st.name)));
    S = load(fullfile(resDir, f(1).name), 'spec', 'T', 'theta1', 'theta2', ...
        'points', 'finalObj', 'baseXi0', 'baseXi1', 'baseXi2', 'baseWraps', ...
        'basePressure');
    spec = S.spec;
    fprintf('\n===== %s (%s) =====\n', st.name, f(1).name);

    % --- regression gate: recorded 5x5 grid, gain 1, 620 kPa --------------
    objRec = BiPam_pulley(spec.name, S.points, spec.cross, 20, S.T, ...
        S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, 0.025, ...
        S.basePressure, S.baseXi0, S.baseXi1, S.baseXi2, S.baseWraps, ...
        struct('nPulleyBPA', 1, 'tackleLineParts', [1 1], ...
        'routingMode', 'moving_via'));
    dmax = 0;
    for k = 1:2
        dmax = max(dmax, max(abs(objRec.Torque_ins(:) - ...
            S.finalObj.Torque_ins(:)), [], 'omitnan'));
    end
    fprintf('REGRESSION recorded-grid rebuild: max |dTorque_ins| = %.3e N*m\n', dmax);
    % off-axis check: |axis projection| vs vector norm at the recorded grid
    for k = 1:2
        axk = spec.kin.axes(k, :) / norm(spec.kin.axes(k, :));
        tv = squeeze(objRec.Torque_ins(:, :, :, k));        % N1x3xN2
        tax = reshape(sum(reshape(tv, [], 3) .* axk, 2), ...
            size(S.T, 3), size(S.T, 4));
        tnv = squeeze(vecnorm(tv, 2, 2));
        rat = abs(tax(:)) ./ tnv(:);
        fprintf('off-axis check crossing %d: min |tau_ax|/|tau_vec| = %.4f (max |tau_ax| %.2f, max |tau_vec| %.2f)\n', ...
            k, min(rat(rat > 0)), max(abs(tax(:))), max(tnv(:)));
    end

    % --- dense grids --------------------------------------------------------
    [Td, th1d, th2d] = buildBiTransforms(spec, N1, N2);
    kinslice = spec.kin;
    kinslice.thetaRangesDeg(2, :) = [0 0];      % secondary pinned at 0 deg
    % NOTE: the BiPamData storage convention stores the two joint slabs on
    % dim 4 (T(:,:,i,1), T(:,:,i,2)), so dim 4 must carry >= 2 pages; the
    % minimal slice grid is N2 = 2 with theta2 = [0 0] (duplicated page,
    % column 1 consumed below).
    [Tm, th1m] = buildBiTransforms(struct('kin', kinslice), N1, 2);
    th1d = th1d(:) * 180 / pi;  th2d = th2d(:) * 180 / pi;  th1m = th1m(:) * 180 / pi;
    % dense-T check against the recorded shared cells
    i_sh = round(linspace(1, N1, 5));
    dT = 0;
    for a = 1:5
        dT = max(dT, max(abs(Td(:, :, i_sh(a), 1) - S.T(:, :, a, 1)), [], 'all'));
    end
    fprintf('dense-grid T check vs recorded shared cells: max |dT| = %.3e\n', dT);

    ax1 = spec.kin.axes(1, :) / norm(spec.kin.axes(1, :));

    % --- required torque on the primary grid --------------------------------
    if st.reqCol > 0
        R = readmatrix(st.reqFile);
        reqAng = R(:, 2);  reqTau = R(:, st.reqCol);
    else
        reqAng = Mmono.ctx.humanAngleD;  reqTau = Mmono.ctx.humanTorque;
    end
    tauReqM = interp1(reqAng, reqTau, th1m, 'linear', 'extrap');

    % --- sweeps ---------------------------------------------------------------
    nG = numel(gainList); nP = numel(pressList);
    tauDense = nan(N1, N2, 2, nG, nP);      % signed hinge-axis torque, dense grid
    infeasD = false(N1, N2, nG, nP);
    slackD  = false(N1, N2, nG, nP);
    tauSlice = nan(N1, nG, nP);             % 25x1 slice, secondary at 0 deg
    infeasS = false(N1, nG, nP);
    slackS  = false(N1, nG, nP);
    buildSec = 0;
    for g = 1:nG
        for pp = 1:nP
            t0 = tic;
            od = BiPam_pulley(spec.name, S.points, spec.cross, 20, Td, ...
                S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, ...
                0.025, pressList(pp), S.baseXi0, S.baseXi1, S.baseXi2, ...
                S.baseWraps, struct('nPulleyBPA', 1, ...
                'tackleLineParts', gainList{g}, 'routingMode', 'moving_via'));
            buildSec = buildSec + toc(t0);
            for k = 1:2
                axk = spec.kin.axes(k, :) / norm(spec.kin.axes(k, :));
                tv = squeeze(od.Torque_ins(:, :, :, k));         % N1x3xN2
                tauDense(:, :, k, g, pp) = reshape( ...
                    sum(reshape(tv, [], 3) .* axk, 2), N1, N2);
            end
            infeasD(:, :, g, pp) = od.PulleyInfeasible;
            slackD(:, :, g, pp) = od.PulleySlack;
            om = BiPam_pulley(spec.name, S.points, spec.cross, 20, Tm, ...
                S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, ...
                0.025, pressList(pp), S.baseXi0, S.baseXi1, S.baseXi2, ...
                S.baseWraps, struct('nPulleyBPA', 1, ...
                'tackleLineParts', gainList{g}, 'routingMode', 'moving_via'));
            tv = squeeze(om.Torque_ins(:, :, 1, 1));              % N1x3
            tauSlice(:, g, pp) = sum(tv .* ax1, 2);
            infeasS(:, g, pp) = om.PulleyInfeasible(:, 1);
            slackS(:, g, pp) = om.PulleySlack(:, 1);
        end
    end
    fprintf('sweeps: %d dense + %d slice builds, %.1f s dense total\n', ...
        nG * nP, nG * nP, buildSec);

    % --- margins + RoM per (G, p) on the primary slice ----------------------
    marginS = nan(N1, nG, nP);
    for g = 1:nG
        for pp = 1:nP
            tr = tauSlice(:, g, pp);
            m = (abs(tr) - abs(tauReqM)) ./ abs(tauReqM);
            m(abs(tauReqM) < 0.5) = NaN;      % guard: near-zero requirement
            m(isnan(tr)) = -Inf;              % infeasible cell: cannot meet
            marginS(:, g, pp) = m;
        end
    end
    for g = 1:nG
        for pp = 1:nP
            pk = max(abs(tauSlice(:, g, pp)), [], 'omitnan');
            mn = min(marginS(:, g, pp), [], 'omitnan');
            rom = romWidth(th1m, marginS(:, g, pp) >= 0);
            [romC, ~] = romContig(th1m, marginS(:, g, pp) >= 0);
            summaryRows{end + 1, 1} = st.name;                  %#ok<AGROW>
            summaryRows{end, 2} = st.joint;
            summaryRows{end, 3} = gainTags{g};
            summaryRows{end, 4} = sprintf('%g', pressList(pp));
            summaryRows{end, 5} = sprintf('%.3f', pk);
            summaryRows{end, 6} = sprintf('%.2f', 100 * mn);
            summaryRows{end, 7} = sprintf('%.2f', rom);
            summaryRows{end, 8} = sprintf('%.2f', romC);
            summaryRows{end, 9} = sprintf('%d', sum(infeasS(:, g, pp)));
            summaryRows{end, 10} = sprintf('%d', sum(slackS(:, g, pp)));
            summaryRows{end, 11} = st.reqSrc;
            fprintf(['  G=%s p=%g kPa: peak|tau|=%7.2f N*m  minMargin=%8.1f%%  ' ...
                'RoM=%6.1f deg (contig %6.1f)  infeas %d/%d\n'], gainTags{g}, ...
                pressList(pp), pk, 100 * mn, rom, romC, ...
                sum(infeasS(:, g, pp)), N1);
        end
    end

    % --- dense pressure axis for the tau(theta, p) surface figure -----------
    % (G = 1/1 slice only; 150 kPa = the valve manifold's regulated floor,
    % 620 kPa = the pressure ceiling of the shared actuator plan)
    pressSurf = [150:50:600, 620];
    tauSurfP = nan(N1, numel(pressSurf));
    for pp = 1:numel(pressSurf)
        os = BiPam_pulley(spec.name, S.points, spec.cross, 20, Tm, ...
            S.finalObj.RestingL, S.finalObj.Kmax, S.finalObj.TendonL, ...
            0.025, pressSurf(pp), S.baseXi0, S.baseXi1, S.baseXi2, ...
            S.baseWraps, struct('nPulleyBPA', 1, ...
            'tackleLineParts', [1 1], 'routingMode', 'moving_via'));
        tv = squeeze(os.Torque_ins(:, :, 1, 1));
        tauSurfP(:, pp) = sum(tv .* ax1, 2);
    end

    % --- save per-site mat ---------------------------------------------------
    out = struct();
    out.site = st.name;
    out.joint = st.joint;
    out.sourceMat = f(1).name;
    out.sizes = struct('RestingL', S.finalObj.RestingL, ...
        'Kmax', S.finalObj.Kmax, 'TendonL', S.finalObj.TendonL);
    out.routingMode = 'moving_via';   % batch mode of record; bowden never used
    out.gainList = gainList;
    out.pressList = pressList;
    out.theta1Dense_deg = th1d;
    out.theta2Dense_deg = th2d;
    out.theta1Slice_deg = th1m;
    out.tauDense = tauDense;
    out.infeasDense = infeasD;
    out.slackDense = slackD;
    out.tauSlice = tauSlice;
    out.infeasSlice = infeasS;
    out.slackSlice = slackS;
    out.tauReqSlice = tauReqM;
    out.reqAngle = reqAng;
    out.reqTau = reqTau;
    out.marginSlice = marginS;
    out.pressSurf = pressSurf;
    out.tauSurfP = tauSurfP;
    out.regressionMaxdT = dmax;
    out.regressionMaxdTorqueIns = dmaxTorque(objRec, S.finalObj);
    out.denseTCheck = dT;
    if strcmp(st.name, 'med_gas_r')
        out.reqMedGasAnkle = reqMedGasAnkle;
        out.medGasAnkleTauDense = tauDense(:, :, 2, :, :);  % crossing 2 = ankle
    end
    if strcmp(st.name, 'rect_fem_r')
        out.routeNote = ['recorded 09-26 route: femur-ending, crossing 2 rigid ' ...
            '(hip live only); gait2392_robotbody.osim since 397a990a reaches ' ...
            'tibia_r with a live knee crossing -- today''s builder yields a ' ...
            'different spec (disclosed, not used here)'];
    end
    matName = fullfile(resDir, sprintf('PulleySurfaces_%s_20261020.mat', st.name));
    save(matName, '-struct', 'out');
    fprintf('saved %s\n', matName);
end

% --- summary CSV -------------------------------------------------------------
csvName = fullfile(resDir, 'PulleySurfaces_summary_20261020.csv');
fid = fopen(csvName, 'w');
hdr = 'site,joint,gain,pressure_kPa,peakAbsTau_Nm,minMargin_pct,romWidth_deg,romContig_deg,nInfeasible_slice,nSlack_slice,reqSource';
fprintf(fid, '%s\n', hdr);
for i = 1:size(summaryRows, 1)
    fprintf(fid, '%s\n', strjoin(string(summaryRows(i, :)), ','));
end
fclose(fid);
fprintf('\nsaved %s\nDONE\n', csvName);

%% =================== local functions =========================================
function [T, theta1, theta2] = buildBiTransforms(spec, N1, N2)
% Verbatim copy of Opt_run_BiPulley.m buildBiTransforms (lines 585-634).
kin = spec.kin;
r1 = rangeOrDefault(kin, 1, [-30, 90]);
r2 = rangeOrDefault(kin, 2, [-110, 20]);
theta1 = linspace(r1(1), r1(2), N1) * pi / 180;
theta2 = linspace(r2(1), r2(2), N2) * pi / 180;

thK = []; xK = []; yK = [];
if isfield(kin, 'kneeSlot') && ~isempty(kin.kneeSlot) && kin.kneeSlot > 0
    if isfield(kin, 'kneeFit') && ~isempty(kin.kneeFit) && ...
            ~isempty(kin.kneeFit.thetaY)
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
% Total angular width (deg) of the set where ok is true: full interval when
% both endpoints ok, half when one (trapezoid measure; NaN-safe because ok
% is a logical with NaN margins mapped to false upstream).
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
% Longest contiguous run where ok is true (deg).
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

function d = dmaxTorque(a, b)
d = 0;
for k = 1:2
    d = max(d, max(abs(a.Torque_ins(:) - b.Torque_ins(:)), [], 'omitnan'));
end
end
