function vasmed_rp3_verify_dense_20261028()
% vasmed_rp3_verify_dense_20261028 -- STUDY (2026-10-28, easteregg2)
%
% Part 3: FRESH-PROCESS DENSE-GRID VERIFICATION of the finalists.
% Coarse 100-pose or 5-node success alone is insufficient (task directive,
% 2026-10-28). This script runs in its own MATLAB process, rebuilds every
% finalist from its saved design vector on a 131-pose near-1-degree grid
% (-119.99994..+9.99999 deg; the -120:1:10 human dense surface is
% interpolated with endpoint-linear extrapolation over the 6e-5 deg
% truncation gap), evaluates through MonoPam_pulley + the X3 fictive-force
% correction, and grades:
%   - worst margin over the ENTIRE dense domain (slack/infeasible/zero-
%     force poses are full shortfalls),
%   - plain-target and 1.05x (route-design 5% margin) meets at the 5
%     surface nodes,
%   - reproducibility: dense worst margin vs the 100-pose search-grid
%     value from rp1/rp2.
%
% Finalists:
%   (0) saved design of record AS-IS (direct drive 1,1, sizes of record)
%   (1) best FIXED-ROUTE sizes-only regrade cell from rp1 (the "before")
%   (2..) every attachment-optimized design from rp2 (the "after" set)
%
% Outputs: vasmed_rp3_dense_20261028.mat, vasmed_rp3_dense_curves.csv,
%          vasmed_rp3_verify_20261028.txt
%
% Invocation: matlab -batch "cd(<folder>); vasmed_rp3_verify_dense_20261028"

R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
here = fileparts(mfilename('fullpath'));
logF = fopen(fullfile(here, 'vasmed_rp3_verify_20261028.txt'), 'w');

C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'ctxCamBest', 'XiUsed');
ctx100 = C.ctxCamBest;
Xi = C.XiUsed(:).';
xRec = C.xFullBest(:).';
p1Rec = xRec(1:3); pEndRec = xRec(4:6);
restRec = xRec(7); tendonRec = xRec(8);
ringOld = ctx100.routeSeed(6, :);
Dia = 20; Pres = 620; KMAX = 0.255; fitting = 0.021; wraps = 3;
FmaxOf = @(rest) maxBPAforce(rest, '20');
cableOff = 0.00075;
cylUp = struct('c', ctx100.geo.tibiaUpperCenter, 'r', ctx100.geo.tibiaUpperRadius + cableOff);
cylLo = struct('c', ctx100.geo.tibiaLowerCenter, 'r', ctx100.geo.tibiaLowerRadius + cableOff);

% ---- dense context: transforms rebuilt from the VERBATIM tables of
% ---- record in buildKneeExtContext20mm.m (copied, not re-derived).
ctxD = ctx100;
[ctxD, phiDense] = denseCtx(ctx100);

TD = readtable(fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022', ...
    'vas_med_surface_dense_20261027.csv'), 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char');
[~, o] = sort(TD.knee_angle_deg(:));
dDeg = TD.knee_angle_deg(o); dTau = TD.tau_h_Nm(o); dArm = TD.arm_m(o);
tauH = interp1(dDeg, dTau, phiDense, 'linear', 'extrap');
armA = interp1(dDeg, dArm, phiDense, 'linear', 'extrap');
T5 = readtable(fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022', ...
    'vas_med_surface_5node_20261027.csv'), 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char');
[~, o5] = sort(T5.knee_angle_deg(:));
nodeDeg = T5.knee_angle_deg(o5); tauHnode = T5.tau_h_Nm(o5);
say(logF, 'Dense verify grid: %d poses %.5f..%.5f deg\n', numel(phiDense), ...
    phiDense(1), phiDense(end));

% ---- assemble finalist list ----
M1 = load(fullfile(here, 'vasmed_rp1_regrade_20261028.mat'), 'rows');
rows1 = M1.rows;
wf = [rows1{:, 7}];
[~, best] = max(wf);
nB1 = rows1{best, 1}; G1 = rows1{best, 2};
fin = {};
fin(end + 1, :) = {'saved_as_is_1_1', 1, 1, [p1Rec, pEndRec, restRec, tendonRec], NaN}; %#ok<AGROW>
fin(end + 1, :) = {sprintf('fixedroute_nB%d_G%d', nB1, G1), nB1, G1, ...
    [p1Rec, pEndRec, rows1{best, 3}, rows1{best, 4}], rows1{best, 5}; %#ok<AGROW>
say(logF, 'rp1 best fixed-route cell: nB%d G%d (search worst margin %+.2f)\n', ...
    nB1, G1, rows1{best, 5});
M2 = load(fullfile(here, 'vasmed_rp2_attachopt_20261028.mat'), 'res');
for k = 1:size(M2.res, 1)
    fin(end + 1, :) = {sprintf('attachopt_nB%d_G%d', M2.res{k, 1}, M2.res{k, 2}), ...
        M2.res{k, 1}, M2.res{k, 2}, M2.res{k, 3}(:).', M2.res{k, 5}}; %#ok<AGROW>
end

out = {}; curves = {};
for k = 1:size(fin, 1)
    x = fin{k, 4}(:).';
    E = evalDense(x, fin{k, 2}, fin{k, 3}, ctxD, Xi, Dia, Pres, KMAX, ...
        fitting, wraps, phiDense, tauH, armA, ringOld, cylUp, cylLo, ...
        nodeDeg, tauHnode, FmaxOf);
    out(end + 1, :) = {fin{k, 1}, fin{k, 2}, fin{k, 3}, x, E.worstMargin, ...
        E.worstMarginPct, E.worstFrac, E.nMeetPlain, E.nMeet5pct, ...
        E.nSlack, E.nInfeas, fin{k, 5}, E}; %#ok<AGROW>
    say(logF, ['VERIFY %-22s nB%d G%d: worst %+.2f N*m (%+.1f%%) | meet %d/%d ' ...
        '(+5%% %d/%d) | slack %d infeas %d | search-grid worst %s\n'], ...
        fin{k, 1}, fin{k, 2}, fin{k, 3}, E.worstMargin, E.worstMarginPct, ...
        E.nMeetPlain, numel(nodeDeg), E.nMeet5pct, E.nSlack, E.nInfeas, ...
        num2str(fin{k, 5}, '%+.2f'));
    curves.(fin{k, 1}) = E;
end

save(fullfile(here, 'vasmed_rp3_dense_20261028.mat'), 'out', 'curves', 'fin', ...
    'phiDense', 'tauH', 'armA', 'nodeDeg', 'tauHnode');
fidC = fopen(fullfile(here, 'vasmed_rp3_dense_curves.csv'), 'w');
fprintf(fidC, 'design,knee_angle_deg,tau_robot_b_Nm,tau_robot_a_Nm,tau_h_Nm,taut,slack,infeas\n');
for k = 1:size(out, 1)
    E = out{k, 13};
    tag = out{k, 1};
    for i = 1:numel(phiDense)
        fprintf(fidC, '%s,%.5f,%.6f,%.6f,%.6f,%d,%d,%d\n', tag, phiDense(i), ...
            E.tauR(i), E.tauR_a(i), E.tauH(i), E.taut(i), E.slack(i), E.infeas(i));
    end
end
fclose(fidC);
say(logF, 'RP3 DONE\n');
fclose(logF);
end

%% ================= LOCAL FUNCTIONS =====================================
function say(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end

%% ------- dense context: VERBATIM transform tables of record -----------
function [ctxD, phiDense] = denseCtx(ctx100)
% Tables copied verbatim from buildKneeExtContext20mm.m lines 23-40.
% Only the pose grid changes: 131 positions over the same RoM instead
% of 100 (near-1-degree spacing).
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533; 0.197344; 0.337395; 0.490178; 1.52146; 2.0944];
knee_x = [-0.0032; 0.00179; 0.00411; 0.0041; 0.00212; -0.001; -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x, knee_x, 'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533; 0.159149; 2.0944];
knee_y = [-0.4226; -0.4082; -0.399; -0.3976; -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y, knee_y, 'cubicspline');
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam = ([23.30 22.22 21.55 21.09 19.91 18.70 17.48 13.82 10.44 7.60 5.52 4.35 4.16 5.01 7.04 10.47]') / 1000;
fcn3 = fit(knee_angle, knee_x_Pam, 'cubicspline');
knee_y_Pam = ([-416.65 -417.03 -417.19 -417.28 -417.41 -417.41 -417.30 -416.28 -414.36 -411.72 -408.62 -405.32 -402.08 -399.16 -396.85 -395.66]') / 1000;
fcn4 = fit(knee_angle, knee_y_Pam, 'cubicspline');
t1_ICR_x = ([29.66 28.54 27.86 27.40 26.23 25.03 23.81 20.03 16.17 12.34 8.67 5.24 2.04 -1.01 -4.1 -7.58]') / 1000;
fcn13 = fit(knee_angle, t1_ICR_x, 'cubicspline');
t1_ICR_y = ([25.97 25.74 25.61 25.53 25.35 25.19 25.03 24.57 24.04 23.39 22.66 21.93 21.32 20.99 21.2 22.33]') / 1000;
fcn14 = fit(knee_angle, t1_ICR_y, 'cubicspline');

kneeMin = -2.0943951; kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, 131);
[~, pos] = min(abs(phi)); phi(pos) = 0;
phiDense = phi * 180 / pi;

N = numel(phi);
T_Pam = zeros(4, 4, N); T_Pam_inv = zeros(4, 4, N);
T_t1_ICR = zeros(4, 4, N); T_ICR_t1 = zeros(4, 4, N);
for i = 1:N
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam = [cos(phi(i)), -sin(phi(i)), 0; sin(phi(i)), cos(phi(i)), 0; 0, 0, 1];
    T_Pam(:, :, i) = RpToTrans(R_Pam, hipToKnee_Pam');
    T_Pam_inv(:, :, i) = T_Pam(:, :, i) \ eye(4);
    t1toICR = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR');
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR');
end
ctxD = ctx100;
ctxD.T_Pam = T_Pam; ctxD.T_Pam_inv = T_Pam_inv;
ctxD.T_ICR_t1 = T_ICR_t1; ctxD.T_t1_ICR = T_t1_ICR;
ctxD.phi = phi(:); ctxD.phiD = phiDense(:);
ctxD.N = N;
ctxD.idxMaxFlex = 1;
end

%% ------- dense evaluator (X3-corrected, both arm methods) -------------
function E = evalDense(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
    phiDense, tauH, armA, ringOld, cylUp, cylLo, nodeDeg, tauHnode, FmaxOf)
p1 = x(1:3); pEnd = x(4:6); rest = x(7); tendon = x(8);
E = struct('ok', false);
[Location, bendM, info] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx);
mL = zeros(numel(phiDense), 1);
for ii = 1:numel(mL)
    tot = 0;
    for j = 1:size(Location, 1) - 1
        pA = Location(j, :, ii);
        pB = Location(j + 1, :, ii);
        if j + 1 == ctx.CrossPoint
            pB = RowVecTrans(ctx.T_Pam(:, :, ii), pB);
        end
        tot = tot + norm(pA - pB);
    end
    mL(ii) = tot;
end
N = numel(mL);
pulCfg = struct('nPulleyBPA', nB, 'tackleLineParts', G, 'gain', G, ...
    'routingMode', 'moving_via');
obj = MonoPam_pulley('rp3', Location, ctx.CrossPoint, Dia, ctx.T_Pam, ...
    rest, rest * (1 - KMAX), tendon, fitting, Pres, Xi(1), Xi(2), Xi(3), ...
    wraps, pulCfg);
sCon = double(obj.sContraction(:));
slack = logical(obj.PulleySlack(:));
infeas = logical(obj.PulleyInfeasible(:));
bend = bendM(:); if numel(bend) ~= N, bend = zeros(N, 1); end
Lm0 = mL - tendon - 2 * fitting - Xi(1);
strain0 = (rest - Lm0) / rest;
comp = max(0, 1 - strain0 / KMAX);
deltaL = Xi(4) .* bend .* comp .^ 2;
Fsingle = max(festo4(Dia, (sCon + deltaL) / rest / KMAX, Pres) .* FmaxOf(rest), 0);
Ftendon = nB * Fsingle / G;
Ftendon(slack) = 0;
Ftendon(infeas) = NaN;
armB = cableArmGrid(pEnd, ringOld, cylUp, cylLo, phiDense);
taut = ~slack & ~infeas & isfinite(Ftendon) & Ftendon > 0;
tauR = Ftendon .* armB; tauR(~taut) = 0;
tauR_a = Ftendon .* armA; tauR_a(~taut) = 0;
meetFrac = zeros(N, 1); meetFrac(taut) = tauR(taut) ./ tauH(taut);
E.tauR = tauR; E.tauR_a = tauR_a; E.tauH = tauH; E.taut = taut;
E.slack = slack; E.infeas = infeas; E.Ftendon = Ftendon; E.armB = armB;
E.worstMargin = min(tauR - tauH);
E.worstMarginPct = 100 * (min(tauR ./ tauH) - 1);
E.worstFrac = min(min(meetFrac), 1);
E.nSlack = nnz(slack); E.nInfeas = nnz(infeas);
tauN = interp1(phiDense, tauR, nodeDeg, 'linear', 'extrap');
tauHn = interp1(phiDense, tauH, nodeDeg, 'linear', 'extrap');
E.nMeetPlain = nnz(tauN >= tauHn - 1e-12);
E.nMeet5pct = nnz(tauN >= 1.05 * tauHn - 1e-12);
E.nodeMargins = tauN - tauHn;
E.ok = true;
end

%% ------- cable arm about the measured moving ICR (study1 verbatim) ----
function armCable = cableArmGrid(pEnd, ringOld, cylUp, cylLo, phiD)
armCable = zeros(numel(phiD), 1);
for i = 1:numel(phiD)
    th = phiD(i) * pi / 180;
    ICR_t1 = [t1icrVal('x', th), t1icrVal('y', th), 0];
    uFirst = cableFirstLeg(pEnd(1:2), ringOld(1:2), {cylUp, cylLo});
    r = pEnd - ICR_t1;
    m = cross(r, [uFirst(1), uFirst(2), 0]);
    armCable(i) = m(3);
end
end

function v = t1icrVal(xy, th)
ka = [0.17, 0.09, 0.03, 0.00, -0.09, -0.17, -0.26, -0.52, -0.79, ...
    -1.05, -1.31, -1.57, -1.83, -2.09, -2.36, -2.62];
if strcmp(xy, 'x')
    tab = [29.66, 28.54, 27.86, 27.40, 26.23, 25.03, 23.81, 20.03, ...
        16.17, 12.34, 8.67, 5.24, 2.04, -1.01, -4.10, -7.58] / 1000;
else
    tab = [25.97, 25.74, 25.61, 25.53, 25.35, 25.19, 25.03, 24.57, ...
        24.04, 23.39, 22.66, 21.93, 21.32, 20.99, 21.20, 22.33] / 1000;
end
v = interp1(ka, tab, th, 'pchip');
end

function [u, nC] = cableFirstLeg(pFrom, pTo, circs)
cur = pFrom(:); dest = pTo(:);
used = false(1, numel(circs));
nC = zeros(1, numel(circs));
u = (dest - cur) / norm(dest - cur);
for it = 1:numel(circs)
    d = dest - cur; L = norm(d); u = d / L;
    hit = 0; hitT = inf;
    for k = 1:numel(circs)
        if used(k), continue; end
        rel = circs{k}.c(:) - cur;
        proj = rel' * u;
        if proj <= 0 || proj >= L, continue; end
        perp = norm(rel - proj * u);
        if perp < circs{k}.r && proj < hitT
            hitT = proj; hit = k;
        end
    end
    if hit == 0, break; end
    c = circs{hit}.c(:); r = circs{hit}.r;
    dd = norm(c - cur);
    if dd <= r, break; end
    base = (c - cur) / dd;
    side = sign((dest(1) - cur(1)) * (c(2) - cur(2)) - ...
        (dest(2) - cur(2)) * (c(1) - cur(1)));
    if side == 0, side = 1; end
    alpha = acos(min(max(r / dd, -1), 1));
    rot = [cos(side * alpha), -sin(side * alpha); ...
        sin(side * alpha), cos(side * alpha)];
    u = rot * base;
    cur = cur + u * sqrt(max(dd^2 - r^2, 0));
    used(hit) = true; nC(hit) = 1;
end
end
