function vasmed_rp1_gate_regrade_20261028()
% vasmed_rp1_gate_regrade_20261028 -- STUDY (2026-10-28, easteregg2)
%
% Part 1 of the reverse-pulley vastus medialis example (Section 4.4 opening).
%
% WHAT THIS SCRIPT IS:
%   (a) Equivalence gates for MonoPam_pulley at nBPA=1, G=1 against
%       MonoPamDataExplicit_balanceX3 (the X3 evaluator family of record):
%         GATE A  Xi3 = 0   : class-vs-class, expect < 1e-6 relative on
%                             single-BPA force and z-torque.
%         GATE B  adopted Xi3: MonoPam_pulley solved state + the X3 fictive
%                             force correction (delta_L = Xi3*bend*comp^2,
%                             the exact Contraction_k formula of the X3
%                             class) vs the X3 class itself. The residual
%                             quantifies the structural approximation: the
%                             pulley class has NO Xi3 input (its closure is
%                             installation-pose rigged), so the fictive
%                             term enters the FORCE evaluation on the
%                             solved equilibrium state, not the closure.
%                             Direction is conservative: fictive loss
%                             LOWERS force => lowers torque => lowers
%                             meet margins.
%   (b) Corrected RE-GRADE of the saved design of record
%       (Vas_Pam_20mm_Result_cam_20261014_0416: p1 [0.040 0.137 0],
%       pEnd of record, rest 0.8277, tendon 0.09249) across
%       nBPA {1..5} x tackle G {1..3} with sizes [rest tendon] re-solved
%       per configuration (patternsearch, study1 protocol) -- the FIXED
%       ROUTE, SIZES-ONLY sweep. This supersedes the 10-27 study1 numbers:
%       grading is now X3-corrected, full-domain (all 100 poses), slack
%       and infeasible poses counted as shortfalls (meetFrac = 0), and
%       the reported minimum margin is over the ENTIRE domain, never over
%       successful nodes only. MATLAB min() silently omits NaN -- the
%       10-27 minMargin was computed exactly that way.
%
% BEN'S RULINGS INCORPORATED (2026-10-28, this session):
%   - The BPA bends in this route (femur-cylinder p2 + condyle/clip p3-p5
%     are BPA-side contacts; only the tibia ring p6-p8 and crosspoint are
%     cable) => Xi3 wrapping loss MUST be incorporated, not just disclosed.
%   - The 0.8277 m resting length of the saved design is TOO LONG; it is
%     reported as the baseline, never as the deliverable.
%   - Grading uses the ACTUAL positive delivered tendon force (Ftendon
%     with PulleySlack/PulleyInfeasible), never reconstructed strain_p or
%     nonzero-torque proxies.
%
% Evaluation conventions (carried from study1_vasmed_canon_20261027.m,
% unchanged where stated): adopted extensor chi pick 32 (XiUsed from the
% cam mat), pressure 620 kPa, KMAX 0.255, wraps 3, fitting 0.021 m,
% moving_via routing, PulleyExitIndex = Cross-1. Arm methods kept
% separate: (a) OpenSim ds/dtheta arm; (b) routed cable arm about the
% measured moving ICR (taut-cable construction to the old proximal ring).
% Sign/domain: knee angle -120..+10 deg, flexion NEGATIVE, human vastus
% medialis extension torque POSITIVE (tau_h > 0 everywhere).
%
% Outputs (this folder):
%   vasmed_rp1_regrade_20261028.mat
%   vasmed_rp1_regrade_summary_20261028.csv
%   vasmed_rp1_gate_20261028.txt  (gate log, also stdout)
%
% Invocation: matlab -batch "cd(<this folder>); vasmed_rp1_gate_regrade_20261028"

R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
here = fileparts(mfilename('fullpath'));
capDir = fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022');
logFile = fullfile(here, 'vasmed_rp1_gate_20261028.txt');
fid = fopen(logFile, 'w');
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

%% ---------------- Fixed data of record --------------------------------
C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'ctxCamBest', 'XiUsed');
xRec = C.xFullBest(:).';           % [p1(1:3) pEnd(1:3) rest tendon]
p1Rec = xRec(1:3); pEndRec = xRec(4:6);
restRec = xRec(7); tendonRec = xRec(8);
ctx = C.ctxCamBest;
Xi = C.XiUsed(:).';                % [Xi0 Xi1 Xi2 Xi3] adopted pick 32
ringOld = ctx.routeSeed(6, :);     % OLD PROXIMAL RING p6 (t1 frame)
Dia = 20; Pres = 620; KMAX = 0.255; fitting = 0.021; wraps = 3;
say(fid, 'Design of record: p1 [%s] pEnd [%s] rest %.6f tendon %.6f\n', ...
    num2str(p1Rec, 5), num2str(pEndRec, 5), restRec, tendonRec);
say(fid, 'Adopted chi (pick 32): Xi0 %.6g Xi1 %.6g Xi2 %.6g Xi3 %.6g\n', Xi);
say(fid, 'Route seed rows 1/9: p1seed [%s] pEndseed [%s]\n', ...
    num2str(ctx.routeSeed(1, :), 5), num2str(ctx.routeSeed(9, :), 5));
say(fid, ['geo.tendonMin (route ctx) = %.4f; this study tendon floor = ' ...
    '0.075 m (pulley packaging)\n'], ctx.geo.tendonMin);

T5 = readtable(fullfile(capDir, 'vas_med_surface_5node_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
[~, ord] = sort(T5.knee_angle_deg(:));
nodeDeg = T5.knee_angle_deg(ord);
tauHnode = T5.tau_h_Nm(ord);
armOsimNode = T5.arm_m(ord);
TD = readtable(fullfile(capDir, 'vas_med_surface_dense_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
[~, ordD] = sort(TD.knee_angle_deg(:));
denseDeg = TD.knee_angle_deg(ordD);          % -120:1:10
tauHdense = TD.tau_h_Nm(ordD);
armOsimDense = TD.arm_m(ordD);
phiD = ctx.phiD(:);                           % 100-pose solver grid
tauHgrid = interp1(denseDeg, tauHdense, phiD, 'linear', 'extrap');
armOsimGrid = interp1(denseDeg, armOsimDense, phiD, 'linear', 'extrap');
say(fid, 'Dense human target: %d poses %g..%g deg, tau_h %.2f..%.2f N*m\n', ...
    numel(denseDeg), denseDeg(1), denseDeg(end), min(tauHdense), max(tauHdense));

%% ---------------- Shared route cache ----------------------------------
routeCache = containers.Map('KeyType', 'double', 'ValueType', 'any');
[Location0, ~, info0] = getRoute(tendonRec, p1Rec, pEndRec, ctx, routeCache); %#ok<ASGLU>
mL0 = routeLengths(Location0, ctx);
Lchar = max(mL0);
say(fid, 'Lchar (max rigid route length) = %.6f m\n', Lchar);

%% ---------------- Arm (b): routed cable about moving ICR ---------------
cableOff = 0.00075;
cylUp = struct('c', ctx.geo.tibiaUpperCenter, 'r', ctx.geo.tibiaUpperRadius + cableOff);
cylLo = struct('c', ctx.geo.tibiaLowerCenter, 'r', ctx.geo.tibiaLowerRadius + cableOff);
armCableRec = cableArmGrid(pEndRec, ringOld, cylUp, cylLo, phiD);

%% ================= GATE A / GATE B ====================================
% Route of record, sizes of record, single BPA, tackle disabled.
pulCfg = struct('nPulleyBPA', 1, 'tackleLineParts', 1, 'gain', 1, ...
    'routingMode', 'moving_via');
obj11 = MonoPam_pulley('gate', Location0, ctx.CrossPoint, Dia, ctx.T_Pam, ...
    restRec, restRec * (1 - KMAX), tendonRec, fitting, Pres, ...
    Xi(1), Xi(2), Xi(3), wraps, pulCfg);

kmaxRec = restRec * (1 - KMAX);
bendRec = info0_bend(routeCache, tendonRec);

% X3 class of record, Xi3 = 0 (GATE A)
x3_0 = MonoPamDataExplicit_balanceX3('gateA', Location0, ctx.CrossPoint, Dia, ...
    ctx.T_Pam, restRec, kmaxRec, tendonRec, fitting, Pres, ...
    Xi(1), Xi(2), Xi(3), 0, wraps, phiD, 1, bendRec);
Fpul = double(obj11.FsingleBPA(:));
Fx3_0 = vecnorm(double(x3_0.F_p), 2, 2);
gateA = max(abs(Fpul - Fx3_0) ./ max(1, abs(Fx3_0)));
tzPul = double(obj11.Torque_p(:, 3));
tzX3_0 = double(x3_0.Torque_p(:, 3));
gateA_t = max(abs(tzPul - tzX3_0) ./ max(1, abs(tzX3_0)));
say(fid, 'GATE A (1,1; Xi3=0): force max rel err %.3g, torque-z max rel err %.3g\n', ...
    gateA, gateA_t);

% X3 class of record, adopted Xi3 (GATE B) vs pulley + post-hoc X3 force
x3_1 = MonoPamDataExplicit_balanceX3('gateB', Location0, ctx.CrossPoint, Dia, ...
    ctx.T_Pam, restRec, kmaxRec, tendonRec, fitting, Pres, ...
    Xi(1), Xi(2), Xi(3), Xi(4), wraps, phiD, 1, bendRec);
Fx3_1 = vecnorm(double(x3_1.F_p), 2, 2);
[~, deltaL] = x3DeltaL(restRec, tendonRec, fitting, mL0, Xi, KMAX, bendRec);
Fcorr = festo4(Dia, (double(obj11.sContraction(:)) + deltaL) / restRec / KMAX, Pres) ...
    .* maxBPAforce(restRec, '20');
Fcorr = max(Fcorr, 0);
gateB = max(abs(Fcorr - Fx3_1) ./ max(1, abs(Fx3_1)));
say(fid, ['GATE B (1,1; adopted Xi3=%.4g): post-hoc X3 force vs X3 class: ' ...
    'max rel err %.3g (mean %.3g); mean delta_L %.4f m, max %.4f m\n'], ...
    Xi(4), gateB, mean(abs(Fcorr - Fx3_1) ./ max(1, abs(Fx3_1))), ...
    mean(deltaL), max(deltaL));

%% ================= Corrected evaluator ================================
% E = evalPulley(x2, nB, G, armMode): x2 = [rest tendon], fixed route of
% record. X3-corrected force on the solved equilibrium state (GATE B
% convention). Full-domain grading on phiD; node grading on nodeDeg.
FmaxOf = @(rest) maxBPAforce(rest, '20');

% Saved-design-as-is regrade (nB=1, G=1, sizes of record) first:
Erec = evalPulley([restRec, tendonRec], 1, 1, 'b', ctx, Xi, Dia, Pres, ...
    KMAX, fitting, wraps, phiD, tauHgrid, armOsimGrid, armCableRec, ...
    p1Rec, pEndRec, ringOld, cylUp, cylLo, routeCache, nodeDeg, ...
    tauHnode, FmaxOf);
printGrade('saved design as-is (1,1)', Erec, tauHnode);

%% ---------------- Sizes-only sweep, 15 configs -------------------------
lb2 = [0.30 * Lchar, 0.075];
ub2 = [0.90 * Lchar, 0.30];
seeds = [0.60 * Lchar, 0.10; 0.45 * Lchar, 0.09; restRec, tendonRec; ...
    0.35, 0.08; 0.55, 0.12];
opts = optimoptions('patternsearch', 'Display', 'none', ...
    'MeshTolerance', 1e-4, 'StepTolerance', 1e-4, 'MaxFunctionEvaluations', 150);
rows = {}; tAll = tic;
for nB = 1:5
    for G = 1:3
        fSeed = inf(size(seeds, 1), 1);
        for sIdx = 1:size(seeds, 1)
            if seeds(sIdx, 2) > ub2(2) || seeds(sIdx, 1) > ub2(1)
                fSeed(sIdx) = 1e6; continue
            end
            E = evalPulley(seeds(sIdx, :), nB, G, 'b', ctx, Xi, Dia, Pres, ...
                KMAX, fitting, wraps, phiD, tauHgrid, armOsimGrid, ...
                armCableRec, p1Rec, pEndRec, ringOld, cylUp, cylLo, ...
                routeCache, nodeDeg, tauHnode, FmaxOf);
            fSeed(sIdx) = -E.worstFrac * 1000 - E.meanFrac * 100;
        end
        [~, ixSeed] = min(fSeed);
        objF = @(x) gradeObj(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, ...
            wraps, phiD, tauHgrid, armOsimGrid, p1Rec, pEndRec, ringOld, ...
            cylUp, cylLo, routeCache, nodeDeg, tauHnode, FmaxOf);
        [xB, ~] = patternsearch(objF, seeds(ixSeed, :), [], [], [], [], ...
            lb2, ub2, [], opts);
        E = evalPulley(xB, nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, fitting, ...
            wraps, phiD, tauHgrid, armOsimGrid, armCableRec, p1Rec, pEndRec, ...
            ringOld, cylUp, cylLo, routeCache, nodeDeg, tauHnode, FmaxOf);
        Ea = evalPulley(xB, nB, G, 'a', ctx, Xi, Dia, Pres, KMAX, fitting, ...
            wraps, phiD, tauHgrid, armOsimGrid, armCableRec, p1Rec, pEndRec, ...
            ringOld, cylUp, cylLo, routeCache, nodeDeg, tauHnode, FmaxOf);
        rows(end + 1, :) = {nB, G, xB(1), xB(2), E.worstMargin, ...
            E.worstMarginPct, E.worstFrac, E.meanFrac, E.meetPlain, ...
            E.meet5pct, E.nMeetPlain, E.nMeet5pct, E.nSlack, E.nInfeas, ...
            E.nTaut, numel(phiD), E.nodeMargins, Ea.worstMargin, ...
            Ea.worstFrac, Ea.nMeetPlain}; %#ok<AGROW>
        say(fid, ['REGRADE nB%d G%d: rest %.4f tendon %.4f | worst margin %+.2f N*m ' ...
            '(%+.1f%%) | meet plain %d/%d, +5%% %d/%d | slack %d infeas %d | ' ...
            'arm(a) worst %+.2f\n'], nB, G, xB(1), xB(2), E.worstMargin, ...
            E.worstMarginPct, E.nMeetPlain, numel(nodeDeg), E.nMeet5pct, ...
            E.nSlack, E.nInfeas, Ea.worstMargin);
    end
end
say(fid, 'Regrade sweep wall time: %.1f min\n', toc(tAll) / 60);

hdr = {'nBPA', 'G', 'rest_m', 'tendon_m', 'worstMargin_Nm', ...
    'worstMarginPct', 'worstFrac', 'meanFrac', 'meetPlainAny', ...
    'meet5pctAny', 'nMeetPlain_5node', 'nMeet5pct_5node', 'nSlack', ...
    'nInfeas', 'nTaut', 'nPoses', 'nodeMarginsNm', 'armA_worstMargin_Nm', ...
    'armA_worstFrac', 'armA_nMeetPlain'};
save(fullfile(here, 'vasmed_rp1_regrade_20261028.mat'), 'rows', 'hdr', ...
    'Erec', 'gateA', 'gateA_t', 'gateB', 'deltaL', 'Xi', 'Lchar', ...
    'p1Rec', 'pEndRec', 'restRec', 'tendonRec', 'ringOld', 'phiD', ...
    'tauHgrid', 'denseDeg', 'tauHdense', 'armOsimDense', 'nodeDeg', ...
    'tauHnode', 'armOsimNode');

fidC = fopen(fullfile(here, 'vasmed_rp1_regrade_summary_20261028.csv'), 'w');
fprintf(fidC, '%s\n', strjoin(hdr(1:end - 5), ','));
for k = 1:size(rows, 1)
    fprintf(fidC, '%d,%d,%.6g,%.6g,%.6g,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%d,%d,%d\n', ...
        rows{k, 1:16});
end
fclose(fidC);
say(fid, 'RP1 DONE\n');
end

%% ================= LOCAL FUNCTIONS =====================================

%% ------- objective wrapper for the sizes-only search -------------------
function f = gradeObj(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
    phiD, tauHgrid, armOsimGrid, p1Rec, pEndRec, ringOld, cylUp, cylLo, ...
    routeCache, nodeDeg, tauHnode, FmaxOf)
E = evalPulley(x, nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
    phiD, tauHgrid, armOsimGrid, [], p1Rec, pEndRec, ringOld, cylUp, cylLo, ...
    routeCache, nodeDeg, tauHnode, FmaxOf);
if ~E.ok
    f = 1e6; return
end
f = -1000 * E.worstFrac - 100 * E.meanFrac;
end

%% ------- the corrected equilibrium evaluator --------------------------
function E = evalPulley(x, nB, G, armMode, ctx, Xi, Dia, Pres, KMAX, ...
    fitting, wraps, phiD, tauHgrid, armOsimGrid, armCable, p1, pEnd, ...
    ringOld, cylUp, cylLo, routeCache, nodeDeg, tauHnode, FmaxOf)
rest = x(1); tendon = x(2);
E = struct('ok', false, 'worstFrac', 0, 'meanFrac', 0, 'worstMargin', -inf, ...
    'worstMarginPct', -100);
if tendon < 0.075 - 1e-12 || tendon > 0.3 || rest <= 0
    return
end
try
    [Location, bendM, info] = getRoute(tendon, p1, pEnd, ctx, routeCache);
catch
    return
end
if isfield(info, 'tendonMax') && tendon > min(info.tendonMax) + 1e-9
    return                       % hard reject (nonlcon cTendonMax convention)
end
mL = routeLengths(Location, ctx);
N = numel(mL);

% --- class of record: equilibrium WITHOUT Xi3 (structural) -------------
pulCfg = struct('nPulleyBPA', nB, 'tackleLineParts', G, 'gain', G, ...
    'routingMode', 'moving_via');
obj = MonoPam_pulley('rp', Location, ctx.CrossPoint, Dia, ctx.T_Pam, ...
    rest, rest * (1 - KMAX), tendon, fitting, Pres, Xi(1), Xi(2), Xi(3), ...
    wraps, pulCfg);
sCon = double(obj.sContraction(:));
slack = logical(obj.PulleySlack(:));
infeas = logical(obj.PulleyInfeasible(:));

% --- X3 fictive-force correction (Ben's ruling 2026-10-28) -------------
bend = bendM(:);
if numel(bend) ~= N, bend = zeros(N, 1); end
[comp, deltaL] = x3DeltaL(rest, tendon, fitting, mL, Xi, KMAX, bend);
Fsingle = festo4(Dia, (sCon + deltaL) / rest / KMAX, Pres) .* FmaxOf(rest);
Fsingle = max(Fsingle, 0);
Ftendon = nB * Fsingle / G;
Ftendon(slack) = 0;
Ftendon(infeas) = NaN;

% --- arm method --------------------------------------------------------
if strcmp(armMode, 'a')
    arm = armOsimGrid;
else
    if isempty(armCable)
        armCable = cableArmGrid(pEnd, ringOld, cylUp, cylLo, phiD);
    end
    arm = armCable;
end

% --- grading: taut + positive delivered force --------------------------
taut = ~slack & ~infeas & isfinite(Ftendon) & Ftendon > 0;
tauR = Ftendon .* arm;
tauR(~taut) = 0;                  % no delivered torque => full shortfall
meetFrac = zeros(N, 1);
meetFrac(taut) = tauR(taut) ./ tauHgrid(taut);
E.tauR = tauR; E.tauH = tauHgrid; E.taut = taut;
E.Ftendon = Ftendon; E.Fsingle = Fsingle; E.deltaL = deltaL; E.comp = comp;
E.worstFrac = min(min(meetFrac), 1);
E.meanFrac = mean(min(meetFrac, 1));
E.worstMargin = min(tauR - tauHgrid);
E.worstMarginPct = 100 * (min(tauR ./ tauHgrid) - 1);
E.nSlack = nnz(slack); E.nInfeas = nnz(infeas); E.nTaut = nnz(taut);
E.slack = slack; E.infeas = infeas;

% node grading (5 surface nodes)
infN = interp1(phiD, double(infeas), nodeDeg, 'linear', 'extrap') > 0.5;
slN = interp1(phiD, double(slack), nodeDeg, 'linear', 'extrap') > 0.5;
tauN = interp1(phiD, tauR, nodeDeg, 'linear', 'extrap');
tauN(infN | slN) = 0;
tauHn = interp1(phiD, tauHgrid, nodeDeg, 'linear', 'extrap');
E.nodeMargins = tauN - tauHn;
meetP = tauN >= tauHn - 1e-12;
meet5 = tauN >= 1.05 * tauHn - 1e-12;
E.nMeetPlain = nnz(meetP); E.nMeet5pct = nnz(meet5);
E.meetPlain = any(meetP); E.meet5pct = any(meet5);
E.tauN = tauN;
E.ok = true;
E.x = x; E.nB = nB; E.G = G; E.rest = rest; E.tendon = tendon;
end

%% ------- X3 fictive loss, EXACT Contraction_k formula of record -------
function [comp, deltaL] = x3DeltaL(rest, tendon, fitting, mL, Xi, KMAX, bend)
% First-pass convention of MonoPamDataExplicit_balanceX3:
% Contraction_k(obj, [], Xi0, [], Xi3): Lmt = rigid MuscleLength,
% gama = 0, comp = max(0, 1 - strain0/KMAX), delta_L = Xi3*bend*comp^2.
Lm0 = mL - tendon - 2 * fitting - Xi(1);
strain0 = (rest - Lm0) / rest;
relstrain0 = strain0 / KMAX;
comp = max(0, 1 - relstrain0);
deltaL = Xi(4) .* bend .* comp .^ 2;
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
% Measured t1-ICR table (fcn13/fcn14), verbatim from
% buildKneeExtContext20mm.m:30-40 / study1_vasmed_canon_20261027.m.
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
% Taut-cable first leg, verbatim from study1_vasmed_canon_20261027.m.
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

%% ------- cached route builder -----------------------------------------
function [Location, bendM, info] = getRoute(tendon, p1, pEnd, ctx, routeCache)
key = round(tendon * 1e9);
if isKey(routeCache, key)
    v = routeCache(key);
    Location = v{1}; bendM = v{2}; info = v{3};
    return
end
[Location, bendM, info] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx);
routeCache(key) = {Location, bendM, info};
end

function bend = info0_bend(routeCache, tendonRec)
v = routeCache(round(tendonRec * 1e9));
bend = v{2}(:);
end

%% ------- muscle lengths (muscleLengthNormal verbatim) -----------------
function mL = routeLengths(Location, ctx)
N = size(Location, 3);
mL = zeros(N, 1);
for ii = 1:N
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
end

%% ------- compact grade printer ----------------------------------------
function printGrade(tag, E, tauHnode)
fprintf(['%s: worst margin %+.2f N*m (%+.1f%%), meet plain %d/%d, +5%% %d/%d, ' ...
    'slack %d, infeas %d, node margins [%s] N*m\n'], tag, E.worstMargin, ...
    E.worstMarginPct, E.nMeetPlain, numel(tauHnode), E.nMeet5pct, ...
    E.nSlack, E.nInfeas, num2str(E.nodeMargins.', '%.2f '));
end

%% ------- local logging helper -----------------------------------------
function say(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end
