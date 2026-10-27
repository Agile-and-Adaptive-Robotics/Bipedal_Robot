%audit_Ext20_nine_20261026.m -- P1-a numerical audit of the extensor 20 mm
% physical data reduction on the NINE-point declared dataset (ExtTest20mm_2
% tab rows; first six = the paper's current six, three later points
% -1/-22/-86.5 deg complete the sweep). READ-ONLY over Ben's sources: the
% xlsx, the cam mat of record, and Knee_Extensor_20mm.m's declared arrays.
% Questions answered numerically:
%   Q1 which adjoint convention reproduces the live torque guard
%      (readserialnumbers2 extensor branch) and Knee_Extensor_20mm.m's own
%      compare section (LC negated, force (cosLC,-sinLC), prf=-0.050)?
%   Q2 what do the flexor-template + sign-flip values (current paper six)
%      equal on all nine points, and how far are they from Q1?
%   Q3 does any candidate measured series track the predicted-at-measured-
%      pressure series (Ben: blue dots sit where yellow should be)?
%   Q4 planar AGENTS.md formula F*sin(LC+0.83 deg)*d as a cross-check.
%   Q5 predictions at measured pressures, at 620 kPa, human targets, and
%      the hybrid classification, exported as the audit CSV.
root = 'D:/GitHub/Bipedal_Robot';
resDir = fullfile(root, 'Code/Matlab/Mesh_Optimization/Results');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
addpath(fullfile(root, 'Code/Matlab/Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data/2022_02_Festo'), '-end');
warning('off', 'all');
cdeg = pi/180;
outdir = fullfile(root, 'tmp/paper_bridge/figscripts');

%% 0. xlsx tab dump (ExtTest20mm_2)
T = readcell(fullfile(root, 'Testing_Data/2026_06_Festo/Results_table_20mm.xlsx'), ...
    'Sheet', 'ExtTest20mm_2');
fprintf('=== ExtTest20mm_2 tab: %d rows x %d cols ===\n', size(T,1), size(T,2));
for r = 1:min(size(T,1), 30)
    rowstr = '';
    for c = 1:min(size(T,2), 22)
        v = T{r,c};
        if isnumeric(v) && isfinite(v)
            rowstr = [rowstr sprintf('%10.4g', v)]; %#ok<AGROW>
        elseif ischar(v)
            rowstr = [rowstr sprintf(' %9.9s', v)]; %#ok<AGROW>
        else
            rowstr = [rowstr '          ']; %#ok<AGROW>
        end
    end
    fprintf('r%02d:%s\n', r, rowstr);
end

%% 1. declared nine-point dataset (Knee_Extensor_20mm.m lines 590-595)
Load9 = [180.2 181.5 168 181 116 100 126 173.4 140];
K9    = [-26.2 -30 -61 -43 -119 5 -1 -22 -86.5];
LC9tab= [28 26 23 27 10.5 21 27 28 22];      % tab entries, positive
P9    = [437 420 425 435 256 480 488 435 350];
n9 = numel(Load9);

%% 2. cam mat of record
C = load(fullfile(resDir, 'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), ...
    'predBest', 'ctxCamBest', 'XiUsed');
ctx = C.ctxCamBest;
pb = C.predBest;
TorqueCam620 = pb.TorqueZ(:);
phiD = ctx.phiD(:);
phi = ctx.phi(:);
fprintf('\nCam ctx T_t1_ICR n=%d; phiD range [%.2f %.2f] deg\n', ...
    size(ctx.T_t1_ICR,3), min(phiD), max(phiD));

%% 3. reduction machinery, both kinematics-table sources
% (a) cam ctx table (what the probe used), (b) the flexor-style measured
% t1-ICR spline table hard-coded in gen_Ext20mm (from Knee_Flexor_Data_20mm)
positions = 100;
knee_angle_tbl = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
t1_ICR_x = ([29.66 28.54 27.86 27.40 26.23 25.03 23.81 20.03 16.17 12.34 8.67 5.24 2.04 -1.01 -4.1 -7.58]')/1000;
t1_ICR_y = ([25.97 25.74 25.61 25.53 25.35 25.19 25.03 24.57 24.04 23.39 22.66 21.93 22.32 21.99 21.2 22.33]')/1000;
fcn13 = fit(knee_angle_tbl, t1_ICR_x, 'cubicspline');
fcn14 = fit(knee_angle_tbl, t1_ICR_y, 'cubicspline');
kneeMin = -2.0943951; kneeMax = 0.17453293;
phiFlx = linspace(kneeMin, kneeMax, positions);
[~, posF] = min(abs(phiFlx)); phiFlx(posF) = 0;
T_t1_ICR_flx = zeros(4,4,positions);
for i = 1:positions
    T_t1_ICR_flx(:,:,i) = RpToTrans(eye(3), ...
        [fcn13(phiFlx(i)), fcn14(phiFlx(i)), 0]');
end

% how close are the two tables?
txA = []; tyA = [];
for i = 1:size(ctx.T_t1_ICR,3)
    txA(end+1,1) = ctx.T_t1_ICR(1,4,i); %#ok<AGROW>
    tyA(end+1,1) = ctx.T_t1_ICR(2,4,i); %#ok<AGROW>
end
fprintf('cam-ctx t1_ICR x range [%.4f %.4f] m, y range [%.4f %.4f] m\n', ...
    min(txA), max(txA), min(tyA), max(tyA));
fprintf('flx-tbl t1_ICR x range [%.4f %.4f] m, y range [%.4f %.4f] m\n', ...
    min(t1_ICR_x), max(t1_ICR_x), min(t1_ICR_y), max(t1_ICR_y));

K_ang = K9 * cdeg;
prfs = {[-0.050 -0.255 0]', [-0.065 -0.255 0]'};
tabnames = {'camCtx', 'flxTbl'};
tabs = {ctx.T_t1_ICR, phi, T_t1_ICR_flx, phiFlx(:)};

variants = {};  % {name, prfIdx, convIdx(1=compare,2=flxflip)}
variants(end+1,:) = {'V1_compare_prf050', 1, 1};
variants(end+1,:) = {'V2_compare_prf065', 2, 1};
variants(end+1,:) = {'V3_flxTemplFlip_prf065', 2, 2};
variants(end+1,:) = {'V4_flxTemplFlip_prf050', 1, 2};

result = zeros(size(variants,1), numel(tabnames), n9);
LC = LC9tab * cdeg;
for ti = 1:numel(tabnames)
    Ticr = tabs{2*ti-1}; ph = tabs{2*ti};
    for v = 1:size(variants,1)
        p_rf = prfs{variants{v,2}};
        T_t1_rf = RpToTrans(eye(3), p_rf);
        Trk = pagemtimes(TransInv(T_t1_rf), Ticr);
        f15 = fit(ph, squeeze(Trk(1,4,:)), 'cubicspline');
        f16 = fit(ph, squeeze(Trk(2,4,:)), 'cubicspline');
        tx = arrayfun(@(a) f15(a), K_ang);
        ty = arrayfun(@(a) f16(a), K_ang);
        if variants{v,3} == 1
            Fx = Load9 .* cos(LC);  Fy = -Load9 .* sin(LC);
        else
            Fx = -Load9 .* cos(LC); Fy = -Load9 .* sin(LC);
        end
        tau = ty .* Fx - tx .* Fy;   % pure-translation adjoint tau_z
        if variants{v,3} == 2
            tau = -tau;
        end
        result(v, ti, :) = tau;
    end
end

% xlsx row anchors (from the dump above)
xlsxTabTorque9 = [44.81 45.59 42.72 45.24 30.15 25.6 31.49 43.12 35.73];
xlsxRow17 = [50.58 32.8 28.8 36.86 49.95 39.71];   % tests 4-9

% V4 planar AGENTS formula: tau = F*sin(LC+0.83deg)*d, d from prf=-0.050
d050 = hypot(-0.050, -0.255);
V5 = Load9 .* sind(LC9tab + 0.83) * d050;
d065 = hypot(-0.065, -0.255);
V6 = Load9 .* sind(LC9tab + 0.83) * d065;

%% 4. predictions
Name = 'Vastus Medialis Proximal Ring 20mm BPA';
rest = pb.rest; tendon = pb.tendon;
kmax = rest * (1 - pb.KMAX);
LocationAdj = pb.Location; bendMeasureAdj = pb.bendMeasure;
CrossPoint = ctx.CrossPoint; Dia = ctx.Dia; fitting = ctx.fitting;
wraps = ctx.wraps; BPAcount = pb.BPAcount; T_Pam = ctx.T_Pam;

predAP = zeros(1, n9);
uP = unique(P9);
for u = 1:numel(uP)
    idx = P9 == uP(u);
    tmpP = MonoPamDataExplicit_balanceX3(Name, LocationAdj, CrossPoint, Dia, ...
        T_Pam, rest, kmax, tendon, fitting, uP(u), ...
        ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3, wraps, phiD, BPAcount, bendMeasureAdj);
    predAP(idx) = interp1(phiD, tmpP.Torque_p(:,3), K9(idx), 'pchip');
end
pred620 = interp1(phiD, TorqueCam620, K9, 'pchip');
human9 = interp1(ctx.humanAngleD(:), ctx.humanTorque(:), K9, 'pchip');

%% 5. reports
xlsxTabTorque6 = [44.8 45.6 42.7 45.2 30.1 25.6];   % tab row 16 (old basis)
taskT6 = [37.23 38.85 37.99 38.31 30.13 22.86];      % at [-26.2 -30 -61 -43 -119 5]

fprintf('\n=== Q1/Q2 nine-point reductions (N*m, tab order) ===\n');
fprintf('angle(deg):        '); fprintf('%9.1f', K9);    fprintf('\n');
fprintf('load(N):           '); fprintf('%9.1f', Load9); fprintf('\n');
fprintf('LC tab (deg):      '); fprintf('%9.1f', LC9tab);fprintf('\n');
fprintf('P(kPa):            '); fprintf('%9.0f', P9);    fprintf('\n');
for ti = 1:2
    for v = 1:4
        fprintf('%s[%s]: ', variants{v,1}, tabnames{ti});
        fprintf('%9.3f', squeeze(result(v, ti, :)));
        fprintf('\n');
    end
end
fprintf('V5_planar_prf050:   '); fprintf('%9.3f', V5); fprintf('\n');
fprintf('V6_planar_prf065:   '); fprintf('%9.3f', V6); fprintf('\n');

fprintf('\n=== anchors ===\n');
fprintf('xlsx row16 (9, tab own):   '); fprintf('%9.2f', xlsxTabTorque9); fprintf('\n');
fprintf('xlsx row17 (tests 4-9):    '); fprintf('%9.2f', xlsxRow17); fprintf('\n');
fprintf('task six    (tab order):   '); fprintf('%9.2f', taskT6); fprintf('\n');
V3six = result(3, 2, 1:6);
V1six = result(1, 2, 1:6);
V1sixA = result(1, 1, 1:6);
V1nine = squeeze(result(1,2,:));
fprintf('V3 flxTbl first six vs task: maxdiff %.3f\n', max(abs(V3six - taskT6)));
fprintf('V1 flxTbl first six vs xlsx16: maxdiff %.3f\n', max(abs(V1six - xlsxTabTorque6)));
fprintf('V1 camCtx first six vs xlsx16: maxdiff %.3f\n', max(abs(V1sixA - xlsxTabTorque6)));
fprintf('V1 flxTbl nine vs xlsx row16: maxdiff %.3f\n', max(abs(V1nine - xlsxTabTorque9)));
fprintf('V2 nine vs xlsx row16: maxdiff %.3f\n', ...
    max(abs(squeeze(result(2,2,:)) - xlsxTabTorque9)));
for v = 1:4
    cand = reshape(result(v,2,:), 1, []);
    m17 = cand(4:9);
    md = max(abs(m17 - xlsxRow17));
    fprintf('V%d flxTbl tests4-9 vs xlsx row17: maxdiff %.3f\n', v, md);
end
md5 = max(abs(V5(4:9) - xlsxRow17));
fprintf('V5 planar tests4-9 vs xlsx row17: maxdiff %.3f\n', md5);
fprintf('V1 vs V3 (nine, flxTbl): maxdiff %.3f\n', ...
    max(abs(squeeze(result(1,2,:)) - squeeze(result(3,2,:)))));
fprintf('V1 camCtx vs V1 flxTbl (nine): maxdiff %.3f\n', ...
    max(abs(squeeze(result(1,1,:)) - squeeze(result(1,2,:)))));

fprintf('\n=== Q3 series tracking (measured candidates vs pred@P) ===\n');
fprintf('pred@P (kPa-specific force model): '); fprintf('%9.2f', predAP); fprintf('\n');
for v = 1:4
    cand = reshape(result(v,2,:), 1, []);
    resd = cand - predAP;
    fprintf('V%d residual vs pred@P: max %7.2f  rms %7.2f\n', v, ...
        max(abs(resd)), sqrt(mean(resd.^2)));
end

fprintf('\n=== Q5 comparison table (adopted V1 flxTbl; all N*m) ===\n');
V1 = reshape(result(1,2,:), 1, []);
fprintf('angle  load_N  LCdeg  P_kPa  meas  pred@P  pred620  human  meas-pred@P  meas-human  meet620\n');
for k = 1:n9
    fprintf('%6.1f %7.1f %6.1f %6.0f %7.2f %7.2f %7.2f %7.2f %10.2f %10.2f %8d\n', ...
        K9(k), Load9(k), LC9tab(k), P9(k), V1(k), predAP(k), pred620(k), ...
        human9(k), V1(k)-predAP(k), V1(k)-human9(k), pred620(k) >= human9(k));
end
fprintf('measured range [%.2f %.2f]; meet-human (meas): %d of %d\n', ...
    min(V1), max(V1), sum(V1 >= human9), n9);

%% 6. CSV export
fid = fopen(fullfile(outdir, 'audit_Ext20_nine_20261026.csv'), 'w');
fprintf(fid, ['series_id,idx,K_deg,Load_N,LC_tab_deg,LC_convention,Pressure_kPa,', ...
    'reaction_x_m,reaction_y_m,reaction_d_m,reaction_ang_deg,measured_torque_Nm,', ...
    'V2_prf065_Nm,V3_flxflip_prf065_Nm,V4_flxflip_prf050_Nm,V5_planar_prf050_Nm,', ...
    'pred_at_measured_P_Nm,pred_620kPa_Nm,human_target_Nm,meas_minus_predP,', ...
    'meas_minus_human,classification\n']);
for k = 1:n9
    fprintf(fid, 'ExtTest20mm_2,%d,%.1f,%.1f,%.1f,tibia-axis(tab positive),%.0f,%.3f,%.3f,%.4f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%s\n', ...
        k, K9(k), Load9(k), LC9tab(k), P9(k), -0.050, -0.255, d050, ...
        atan2d(-0.255,-0.050), V1(k), result(2,2,k), result(3,2,k), ...
        result(4,2,k), V5(k), predAP(k), pred620(k), human9(k), ...
        V1(k)-predAP(k), V1(k)-human9(k), 'measured(adjoint,guard convention)');
end
fclose(fid);
fprintf('WROTE audit_Ext20_nine_20261026.csv\n');
disp('DONE_audit_Ext20_nine');
