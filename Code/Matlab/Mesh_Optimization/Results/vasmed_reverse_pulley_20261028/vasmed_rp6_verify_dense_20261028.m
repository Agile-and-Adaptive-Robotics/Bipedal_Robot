function vasmed_rp6_verify_dense_20261028()
% vasmed_rp6_verify_dense_20261028 -- FRESH-PROCESS DENSE VERIFICATION
% (2026-10-28). Coarse 100-pose or 5-node success alone is insufficient;
% this script re-evaluates the rp4 finalists from their saved design
% vectors on a 131-pose near-1-degree grid, through the same verified
% closure (BiPam_pulley2 equation per pose, Xi3 inside the equilibrium,
% cb from the verbatim LokBi2 projection), and grades full-domain worst
% margins + node meets (plain and 1.05x) + slack/infeasible counts, and
% reproducibility vs the 100-pose search values.
%
% Finalists: saved design as-is (1,1); best fixed-route sizes-only cell
% from rp4 Stage B; every attachment-optimized design from rp4 Stage C.
%
% Outputs: vasmed_rp6_dense_20261028.mat, vasmed_rp6_dense_curves.csv,
%          vasmed_rp6_verify_20261028.txt
%
% Invocation: matlab -batch "cd(<folder>); vasmed_rp6_verify_dense_20261028"

R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
here = fileparts(mfilename('fullpath'));
logF = fopen(fullfile(here, 'vasmed_rp6_verify_20261028.txt'), 'w');

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

ctxD = denseCtx(ctx100);
phiDense = ctxD.phiD(:);

capDir = fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022');
TD = readtable(fullfile(capDir, 'vas_med_surface_dense_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
[~, oD] = sort(TD.knee_angle_deg(:));
tauH = interp1(TD.knee_angle_deg(oD), TD.tau_h_Nm(oD), phiDense, 'linear', 'extrap');
armA = interp1(TD.knee_angle_deg(oD), TD.arm_m(oD), phiDense, 'linear', 'extrap');
T5 = readtable(fullfile(capDir, 'vas_med_surface_5node_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
[~, o5] = sort(T5.knee_angle_deg(:));
nodeDeg = T5.knee_angle_deg(o5); tauHnode = T5.tau_h_Nm(o5);
say(logF, 'Dense grid: %d poses %.5f..%.5f deg\n', numel(phiDense), ...
    phiDense(1), phiDense(end));

M4 = load(fullfile(here, 'vasmed_rp4_final_20261028.mat'), 'rows', 'res');
rows = M4.rows;
wf = [rows{:, 7}];
[~, best] = max(wf);
fin = {};
fin(end + 1, :) = {'saved_as_is_1_1', 1, 1, [p1Rec, pEndRec, restRec, tendonRec], NaN}; %#ok<AGROW>
fin(end + 1, :) = {sprintf('fixedroute_nB%d_G%d', rows{best, 1}, rows{best, 2}), ...
    rows{best, 1}, rows{best, 2}, [p1Rec, pEndRec, rows{best, 3}, rows{best, 4}], ...
    rows{best, 5}}; %#ok<AGROW>
say(logF, 'rp4 best fixed-route cell: nB%d G%d (search worst margin %+.2f)\n', ...
    rows{best, 1}, rows{best, 2}, rows{best, 5});
for k = 1:size(M4.res, 1)
    fin(end + 1, :) = {sprintf('attachopt_nB%d_G%d', M4.res{k, 1}, M4.res{k, 2}), ...
        M4.res{k, 1}, M4.res{k, 2}, M4.res{k, 3}(:).', M4.res{k, 5}}; %#ok<AGROW>
end

out = {}; curves = struct();
for k = 1:size(fin, 1)
    x = fin{k, 4}(:).';
    E = evalDense(x, fin{k, 2}, fin{k, 3}, ctxD, Xi, Dia, Pres, KMAX, ...
        fitting, wraps, phiDense, tauH, armA, ringOld, cylUp, cylLo, ...
        nodeDeg, FmaxOf);
    out(end + 1, :) = {fin{k, 1}, fin{k, 2}, fin{k, 3}, x, E.worstMargin, ...
        E.worstMarginPct, E.worstFrac, E.nMeetPlain, E.nMeet5pct, ...
        E.nSlack, E.nInfeas, fin{k, 5}, E}; %#ok<AGROW>
    say(logF, ['VERIFY %-22s nB%d G%d: worst %+.2f N*m (%+.1f%%) meet %d/5 ' ...
        '(+5%% %d/5) slack %d infeas %d | 100-pose worst %s\n'], fin{k, 1}, ...
        fin{k, 2}, fin{k, 3}, E.worstMargin, E.worstMarginPct, ...
        E.nMeetPlain, E.nMeet5pct, E.nSlack, E.nInfeas, num2str(fin{k, 5}, '%+.2f'));
    curves.(fin{k, 1}) = E;
end

save(fullfile(here, 'vasmed_rp6_dense_20261028.mat'), 'out', 'curves', 'fin', ...
    'phiDense', 'tauH', 'armA', 'nodeDeg', 'tauHnode');
fidC = fopen(fullfile(here, 'vasmed_rp6_dense_curves.csv'), 'w');
fprintf(fidC, 'design,knee_angle_deg,tau_robot_b_Nm,tau_robot_a_Nm,tau_h_Nm,taut,slack,infeas\n');
for k = 1:size(out, 1)
    E = out{k, 13}; tag = out{k, 1};
    for i = 1:numel(phiDense)
        fprintf(fidC, '%s,%.5f,%.6f,%.6f,%.6f,%d,%d,%d\n', tag, phiDense(i), ...
            E.tauR(i), E.tauR_a(i), E.tauH(i), E.taut(i), E.slack(i), E.infeas(i));
    end
end
fclose(fidC);
say(logF, 'RP6 DONE\n');
fclose(logF);
end

%% ================= LOCAL FUNCTIONS =====================================
function say(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end

function [ctxD, phiDense] = denseCtx(ctx100)
% VERBATIM transform tables from buildKneeExtContext20mm.m lines 23-40;
% only the pose count changes (131 -> near-1-degree spacing).
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
T_t1_ICR = zeros(4, 4, N); T_ICR_t1 = zeros(4, 4, N); T_f_t1 = zeros(4, 4, N);
for i = 1:N
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam = [cos(phi(i)), -sin(phi(i)), 0; sin(phi(i)), cos(phi(i)), 0; 0, 0, 1];
    T_Pam(:, :, i) = RpToTrans(R_Pam, hipToKnee_Pam');
    T_Pam_inv(:, :, i) = T_Pam(:, :, i) \ eye(4);
    t1toICR = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR');
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR');
    T_f_t1(:, :, i) = (T_t1_ICR(:, :, i) * T_Pam_inv(:, :, i)) \ eye(4);
end
ctxD = ctx100;
ctxD.T_Pam = T_Pam; ctxD.T_Pam_inv = T_Pam_inv;
ctxD.T_ICR_t1 = T_ICR_t1; ctxD.T_t1_ICR = T_t1_ICR; ctxD.T_f_t1 = T_f_t1;
ctxD.phi = phi(:); ctxD.phiD = phiDense(:);
ctxD.N = N; ctxD.idxMaxFlex = 1;
end

function E = evalDense(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
    phiDense, tauH, armA, ringOld, cylUp, cylLo, nodeDeg, FmaxOf)
p1 = x(1:3); pEnd = x(4:6); rest = x(7); tendon = x(8);
E = struct('ok', false);
[Location, bend, info] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx);
N = numel(phiDense);
mL = zeros(N, 1); span = zeros(N, 1);
for ii = 1:N
    tot = 0;
    for j = 1:size(Location, 1) - 1
        pA = Location(j, :, ii); pB = Location(j + 1, :, ii);
        if j + 1 == ctx.CrossPoint
            pB = RowVecTrans(ctx.T_Pam(:, :, ii), pB);
        end
        tot = tot + norm(pA - pB);
    end
    mL(ii) = tot;
    raw = info.raw(:, :, ii);
    r5t1 = RowVecTrans(ctx.T_f_t1(:, :, ii), raw(5, :));
    s2 = norm(r5t1 - raw(6, :));
    for r = 6:8
        s2 = s2 + norm(raw(r + 1, :) - raw(r, :));
    end
    span(ii) = s2;
end
bend = bend(:); if numel(bend) ~= N, bend = zeros(N, 1); end
Lm0 = mL - tendon - 2 * fitting - Xi(1);
comp = max(0, 1 - ((rest - Lm0) / rest) / KMAX);
deltaL = Xi(4) .* bend .* comp .^ 2;
s0x3 = rest - mL + Xi(1) + tendon + 2 * fitting + deltaL;
cb = bracketCb(Location, ctx, Xi, N);
P = psiSolve(rest, nB, G, s0x3, span, cb, KMAX, Dia, Pres, tendon, wraps, FmaxOf);
Ftendon = P.Ftendon; slack = P.slack; infeas = P.infeas;
armB = cableArmGrid(pEnd, ringOld, cylUp, cylLo, phiDense);
taut = ~slack & ~infeas & isfinite(Ftendon) & Ftendon > 0;
tauR = Ftendon .* armB; tauR(~taut) = 0;
tauR_a = Ftendon .* armA; tauR_a(~taut) = 0;
meetFrac = zeros(N, 1); meetFrac(taut) = tauR(taut) ./ tauH(taut);
E.tauR = tauR; E.tauR_a = tauR_a; E.tauH = tauH; E.taut = taut;
E.slack = slack; E.infeas = infeas; E.Ftendon = Ftendon;
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

function P = psiSolve(rest, nB, G, s0x3, span, cb, KMAX, Dia, Pres, ...
    tendon, wraps, FmaxOf)
N = numel(s0x3);
[~, iRef] = max(span);
DL = span(iRef) - span;
Fmax = FmaxOf(rest);
cSpr = 1 / (wraps * 1.51e-6 * 193e9 / tendon);
Ffun = @(s) festo4(Dia, (s / rest) / KMAX, Pres) * Fmax;
sHi = KMAX * rest;
sCon = zeros(N, 1); Fsingle = zeros(N, 1);
slack = false(N, 1); infeas = false(N, 1);
for i = 1:N
    psi = @(s) G * s - G * cb(i) * Ffun(s) - nB * Ffun(s) * cSpr / G ...
        - DL(i) - G * (s0x3(i) - DL(i));
    pl = psi(0); ph = psi(sHi);
    if pl > 0
        sCon(i) = 0; Fsingle(i) = Ffun(0); slack(i) = true;
    elseif ph < 0
        sCon(i) = sHi; Fsingle(i) = Ffun(sHi); infeas(i) = true;
    else
        try
            sCon(i) = max(fzero(psi, [0, sHi]), 0);
        catch
            sCon(i) = 0; slack(i) = true;
        end
        Fsingle(i) = Ffun(sCon(i));
    end
end
Ftendon = nB * Fsingle / G;
Ftendon(slack) = 0;
Ftendon(infeas) = NaN;
P = struct('Ftendon', Ftendon, 'Fsingle', Fsingle, 'sCon', sCon, ...
    'slack', slack, 'infeas', infeas, 'DL', DL);
end

function cb = bracketCb(Location, ctx, Xi, N)
cb = zeros(N, 1);
C = diag(1 ./ [Xi(2), Xi(3), Xi(3)]);
pA = Location(1, :, 1);
Pbr = [9.48, -36.21, 30.86] / 1000;
phbrA = pA - Pbr;
thbrA = atan2(phbrA(2), phbrA(1));
RhbrZ = [cos(thbrA), -sin(thbrA), 0; sin(thbrA), cos(thbrA), 0; 0, 0, 1];
pbrhA = RhbrZ' * phbrA';
thetaY = atan2(pbrhA(3), pbrhA(1));
Ry = [cos(thetaY), 0, sin(thetaY); 0, 1, 0; -sin(thetaY), 0, cos(thetaY)];
Thbr = RpToTrans(RhbrZ * Ry', Pbr');
Tinv = Thbr \ eye(4);
for ii = 1:N
    pt1 = Location(1, :, ii);
    pt2 = pt1;
    for r = 2:size(Location, 1)
        if norm(Location(r, :, ii) - pt1) > 1e-6
            pt2 = Location(r, :, ii);
            if r == ctx.CrossPoint
                pt2 = RowVecTrans(ctx.T_Pam(:, :, ii), pt2);
            end
            break
        end
    end
    u = pt2 - pt1;
    if norm(u) < 1e-9
        cb(ii) = 1 / Xi(2);
        continue
    end
    ub = RowVecTrans(Tinv, u / norm(u));
    cb(ii) = (ub * C) * ub.';
end
end

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
