%% sweep_flx20_tendon.m
% Tendon-length sweep for the flexor 2x20 mm BPA pair under KMAX = 0.235
% (Ben, 2026-10-22: the 20 mm BPA max contraction is 23.5%, not 25.5%).
% Question: must the tendon be lengthened or shortened so the pair still
% meets/exceeds the human Bifemsh (short head) torque target at <= 620 kPa?
%
% The kinematics, fits, contexts, and MonoPam_mult calls below are copied
% from Knee_Flexor_Data_20mm.m so the numbers match that script.
%
% Two routing conventions are swept:
%   convB (script torque sections): route built once with xBest(8),
%         TendonL = xBest(8) - 0.02 + dT
%   convA (physical tendon swap):   route and TendonL both use
%         xBest(8) - 0.02 + dT (same pattern as the script's bone-plot
%         section, line ~791)
clc; clear; close all;

here = fileparts(mfilename('fullpath'));
cd(here);
root = here;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', here)
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

%% Joint rotation transformation matrices (verbatim from the script)
positions = 100;

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);
R_Pam = zeros(3, 3, positions);
T_Pam = zeros(4, 4, positions);
t1toICR = zeros(1,3,positions);
T_t1_ICR = zeros(4, 4, positions);
T_ICR_t1 = zeros(4, 4, positions);

c = pi/180;

knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,t1_ICR_y,'cubicspline');

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, positions);
[val, pos] = min(abs(phi));
phi(pos) = 0;

for i = 1:positions
    hipToKnee = [fcn1(phi(i)), fcn2(phi(i)), 0];
    R(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    T(:, :, i) = RpToTrans(R(:, :, i), hipToKnee');

    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');

    t1toICR(1,:,i) = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR(1,:,i)');
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR(1,:,i)');
end

phiD = phi*180/pi;

%% Optimized design from the result mat (same loads as the script)
resultData = load('Bifemsh_20mm_Result.mat', 'xBest', 'routeCtx', 'Xi3');
xBest = resultData.xBest;
routeCtx = resultData.routeCtx;
routeCtx.geo = buildGeoExclusion();

try
    resultXi = load('Bifemsh_20mm_Result.mat', 'XiUsed');
    Xi = resultXi.XiUsed;
    fprintf('XiUsed from mat: %.6g %.6g %.6g %.6g\n', Xi);
catch
    L = load('minimizeFlxPin10_results_20260730_2transforms_Z2.mat', ...
        'filtered_results', 'xCols');
    g = filtered_results(1, xCols);
    Xi = [g(1) g(2) g(3) resultData.Xi3];
    fprintf('Xi from legacy 20260730 pick 1: %.6g %.6g %.6g %.6g\n', Xi);
end
Xi0 = Xi(1); Xi1 = Xi(2); Xi2 = Xi(3); Xi3 = Xi(4);

p1 = xBest(1:3);
p2 = xBest(4:6);
rest = xBest(7);
tendonBase = xBest(8) - 0.02;
fprintf('\nrest (xBest(7)) = %.6f m, xBest(8) = %.6f m, tendon = %.6f m\n', ...
    rest, xBest(8), tendonBase);

Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2;
Dia = 20;
fitting = 0.021;
wraps = 6;
BPAcount = 2;
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);

%% Human Bifemsh target interpolated on the phi grid
H = readmatrix(fullfile(root, 'Testing_Data', '2022_02_Festo', ...
    'OpenSim_Bifem_Results.txt'), 'FileType', 'text', 'NumHeaderLines', 7);
humanAngle = H(:,2);
TorqueHz = H(:,4);
[humanAngleU, iuH] = unique(humanAngle);
humanTgtPhi = interp1(humanAngleU, TorqueHz(iuH), phiD, 'pchip');
humanTgtPhi = humanTgtPhi(:);
demand = isfinite(humanTgtPhi) & humanTgtPhi < 0;
fprintf('angles with flexion demand: %d of %d\n', sum(demand), positions);

%% Routes: convB fixed (xBest(8)) and per-tendon convA built inside the loop
[LocationB_fix, bmA_fix] = buildKneeFlexorRoute20mm(p1, p2, xBest(8), routeCtx);
[LocationB_fix2, bmB_fix] = buildKneeFlexorRoute20mm(p1B, p2B, xBest(8), routeCtx);
LPfix = {LocationB_fix; LocationB_fix2};
BMfix = {bmA_fix; bmB_fix};

%% Sweep
KMAXset = [0.255, 0.235];
dTmm = [-100 -95 -90 -85 -80 -75 -70 -65 -60 -55 -50 -45 -40 -35 -30 ...
        -25 -20 -15 -10 -5 0 5 10 15 20 25 30 35 40];
convNames = {'convA_physical', 'convB_script'};
sweep = [];

for iK = 1:numel(KMAXset)
    KMAX = KMAXset(iK);
    kmax = (1 - KMAX)*rest;
    for iC = 1:2
        for iT = 1:numel(dTmm)
            tendon_s = tendonBase + dTmm(iT)/1000;
            if iC == 1
                [Loc1, bm1] = buildKneeFlexorRoute20mm(p1, p2, tendon_s, routeCtx);
                [Loc2, bm2] = buildKneeFlexorRoute20mm(p1B, p2B, tendon_s, routeCtx);
                LP = {Loc1; Loc2}; BM = {bm1; bm2};
            else
                LP = LPfix; BM = BMfix;
            end
            bpa = MonoPam_mult(Name, LP, CrossPoint, Dia, T_Pam, rest, ...
                kmax, tendon_s, fitting, 620, Xi0, Xi1, Xi2, Xi3, wraps, ...
                phiD, BPAcount, BM);
            tq = bpa.Torque_p(:,3);
            marg = humanTgtPhi(demand) - tq(demand);   % >=0 meets/exceeds
            marg(isnan(marg)) = -inf;
            nFail = sum(marg < 0);
            [minM, im] = min(marg);
            sAll = [bpa.strain_p{1}; bpa.strain_p{2}];
            sweep(end+1,:) = [iK, iC, dTmm(iT), tendon_s*1000, nFail, ...
                minM, phiD(demand(im)), min(sAll(:)), max(sAll(:))]; %#ok<SAGROW>
            fprintf(['KMAX %.3f %-15s dT %+4.0f mm  tendon %6.1f mm  ' ...
                'nFail %3d  minMargin %+8.2f Nm @ %7.1f deg  ' ...
                'strain [%+.4f, %+.4f]\n'], KMAX, convNames{iC}, ...
                dTmm(iT), tendon_s*1000, nFail, minM, phiD(demand(im)), ...
                min(sAll(:)), max(sAll(:)));
        end
    end
end

SweepTable = array2table(sweep, 'VariableNames', {'KMAX_idx', 'conv_idx', ...
    'dTendon_mm', 'tendon_mm', 'nFailAt620kPa', 'minMargin_Nm', ...
    'worstAngle_deg', 'minStrain', 'maxStrain'});
writetable(SweepTable, 'flx20_kmax235_tendon_sweep.csv');
fprintf('\nSweep table: flx20_kmax235_tendon_sweep.csv\n');

%% Pass windows summary
for iK = 1:numel(KMAXset)
    for iC = 1:2
        sel = sweep(sweep(:,1)==iK & sweep(:,2)==iC, :);
        pass = sel(sel(:,5)==0, :);
        if isempty(pass)
            fprintf('KMAX %.3f %-15s: NO passing tendon in [-100,+40] mm\n', ...
                KMAXset(iK), convNames{iC});
        else
            fprintf(['KMAX %.3f %-15s: pass dTendon %+.0f .. %+.0f mm ' ...
                '(tendon %.1f .. %.1f mm)\n'], KMAXset(iK), convNames{iC}, ...
                min(pass(:,3)), max(pass(:,3)), min(pass(:,4)), max(pass(:,4)));
        end
    end
end

%% Pressure back-calculation at the smallest passing shortening, convA
% (closest to the current hardware) and convB, KMAX = 0.235.
Pgrid = (0:5:620)';
for iC = 1:2
    sel = sweep(sweep(:,1)==2 & sweep(:,2)==iC, :);
    pass = sel(sel(:,5)==0, :);
    if isempty(pass), continue; end
    [~, iWin] = max(pass(:,3));           % largest dT (least change) that passes
    rowWin = pass(iWin, :);
    tendon_s = rowWin(4)/1000;
    KMAX = 0.235; kmax = (1-KMAX)*rest;
    TqVsP = nan(numel(Pgrid), positions);
    for kP = 1:numel(Pgrid)
        if iC == 1
            [Loc1, bm1] = buildKneeFlexorRoute20mm(p1, p2, tendon_s, routeCtx);
            [Loc2, bm2] = buildKneeFlexorRoute20mm(p1B, p2B, tendon_s, routeCtx);
            LP = {Loc1; Loc2}; BM = {bm1; bm2};
        else
            LP = LPfix; BM = BMfix;
        end
        tmpP = MonoPam_mult(Name, LP, CrossPoint, Dia, T_Pam, rest, kmax, ...
            tendon_s, fitting, Pgrid(kP), Xi0, Xi1, Xi2, Xi3, wraps, phiD, ...
            BPAcount, BM);
        TqVsP(kP,:) = tmpP.Torque_p(:,3).';
    end
    Preq = nan(positions,1);
    Flag = repmat("no human data", positions, 1);
    for i = 1:positions
        tgt = humanTgtPhi(i);
        if ~isfinite(tgt), continue; end
        if tgt >= 0, Flag(i) = "no flexion demand"; continue; end
        col = TqVsP(:,i);
        fin = isfinite(col);
        if ~any(fin), Flag(i) = "BPA out of range"; continue; end
        if min(col(fin)) > tgt, Flag(i) = "unreachable at 620 kPa"; continue; end
        kk = find(col <= tgt, 1);
        if kk <= 1
            Preq(i) = 0;
        else
            pq = Pgrid([kk-1 kk]);
            [xs, isrt] = sort(col([kk-1 kk]));
            Preq(i) = interp1(xs, pq(isrt), tgt);
        end
        Flag(i) = "ok";
    end
    fprintf(['\nPressure back-calc, KMAX=0.235, %s, dTendon %+.0f mm ' ...
        '(tendon %.1f mm): max required = %.1f kPa\n'], convNames{iC}, ...
        rowWin(3), tendon_s*1000, max(Preq));
    tab = table(phiD(:), humanTgtPhi(:), Preq(:), Flag, ...
        'VariableNames', {'knee_angle_deg','human_target_Nm', ...
        'required_pressure_kPa','status'});
    writetable(tab, sprintf('flx20_kmax235_%s_dT%+.0fmm_pressure.csv', ...
        convNames{iC}, rowWin(3)));
end
fprintf('\nSWEEP DONE\n');
