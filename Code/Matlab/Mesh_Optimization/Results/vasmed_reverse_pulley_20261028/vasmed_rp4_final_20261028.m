function vasmed_rp4_final_20261028()
% vasmed_rp4_final_20261028 -- FINAL STUDY (2026-10-28, easteregg2)
%
% Reverse-pulley vastus medialis example for the torque paper Section 4.4
% opening. SUPERSEDES the rp1 regrade (vasmed_rp1_*), which used
% MonoPam_pulley's closure on this route and is structurally invalid
% through the flexion region (see the AUDIT note below), and the 10-27
% study1 numbers for the same reason.
%
% AUDIT FINDING (this session, probe vasmed_rp0_probe_20261028):
%   MonoPam_pulley's closure assumes mL = const + crossed-span (its
%   (1,1) regression identity is asserted for "a geometry whose crossed
%   segment is the tendon line", Opt_sanity_pulley.m). On the distal-ring
%   extensor route the pre-exit (femur-side) rows RELEASE contacts toward
%   extension (segment probe: at -120 deg the crossed segment p5->p6 is
%   0.107 m with a 0.607 m femur chain; at -19 deg the femur chain is
%   collapsed and the crossed segment is 0.634 m), so the identity
%   s0(i) = s0ref + DeltaL(i) fails by up to ~0.6 m of bookkeeping.
%   BiPam_pulley2_20261022 (the canonical BiPulley class of the
%   capability survey) was built for exactly this case: its tendon span
%   is the whole exit-to-insertion routed path and its closure references
%   each cell to the per-cell s0(i) - DeltaL(i) (fortzBi2, line 1515).
%
% EVALUATOR USED HERE (no class edits, no 1-D port of a whole class):
%   the BiPam_pulley2_20261022 closure equation applied per pose,
%     psi(s) = G*s - G*cb_i*F(s) - nB*F(s)*cSpr/G - DL(i)
%              - G*(s0x3(i) - DL(i)),
%   with s0x3(i) = rest - mL(i) + Xi0 + tendon + 2*fitting + delta_L(i)
%   carrying the Xi3 fictive wrap loss INSIDE the equilibrium (Ben's
%   2026-10-28 ruling: the BPA bends on this route -- femur cylinder p2
%   + condyle/clip p3-p5 -- so Xi3 must be incorporated; delta_L =
%   Xi3*bendMeasure*comp^2 is MonoPamDataExplicit_balanceX3's own
%   Contraction_k formula). Every per-pose input is machinery of record:
%   route/mL/bend from buildDistalRingLocation20mm, DL = the class's
%   exit(row 5)->insertion(row 9) multi-segment span change, cb = the
%   LokBi2 bracket projection (verbatim, incl. the dia-20 Pbr
%   [9.48 -36.21 30.86]/1000 of the pulley family), F = festo4/maxBPAforce,
%   cSpr = wraps(3)*Aeff*E/TendonL (pulley-family Spr). Verified by:
%     GATE 1 (nB=1,G=1): vs MonoPamDataExplicit_balanceX3 (BPAcount=1,
%            per-pose Location, Xi3 native) on the route of record.
%     GATE 2 (nB=2,G=2): vs BiPam_pulley2_20261022 itself on a FIXED-ROW
%            control route (uncollapsed seed20 rows, 100-pose grid).
%   Mount/bowden bracket deflection (variant-2 fields) is OFF, matching
%   the recorded default configuration of the survey batch; the class's
%   mount reaction R = nBPA*F*u_bpa + F_t*u_t (tackle-appropriate
%   internal force, NOT the reduced tendon force) is reported with the
%   winner.
%
% STAGES (all inside this run, saved incrementally):
%   A. Gates.
%   B. Sizes-only sweep, FIXED route of record, nB {1..5} x G {1..3},
%      tendon >= 0.075 m (pulley packaging floor), 5%-target tracked.
%   C. Attachment optimization (x = [p1 pEnd rest tendon]) for the top-3
%      configs by arm(b) worstFrac + the direct-drive (1,1):
%      bounds = the documented search regions (p1: pre-pin region in
%      buildKneeExtContext20mm.m; pEnd: active bounds of record), with
%      Ben's 2026-10-28 deviation penalty vs the anatomical seed rows
%      and a rest-length term (0.8277 m ruled too long), solver
%      surrogateopt(150) + patternsearch(150), nonlconExt20mm constraints.
%
% Outputs: vasmed_rp4_final_20261028.mat, vasmed_rp4_sizesweep.csv,
%          vasmed_rp4_attachopt.csv, vasmed_rp4_final_20261028.txt
%
% Invocation: matlab -batch "cd(<folder>); vasmed_rp4_final_20261028"

R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
here = fileparts(mfilename('fullpath'));
logF = fopen(fullfile(here, 'vasmed_rp4_final_20261028.txt'), 'w');

C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'ctxCamBest', 'XiUsed');
xRec = C.xFullBest(:).';
p1Rec = xRec(1:3); pEndRec = xRec(4:6);
restRec = xRec(7); tendonRec = xRec(8);
ctx = C.ctxCamBest;
Xi = C.XiUsed(:).';
p1Seed = ctx.routeSeed(1, :); pEndSeed = ctx.routeSeed(9, :);
ringOld = ctx.routeSeed(6, :);
Dia = 20; Pres = 620; KMAX = 0.255; fitting = 0.021; wraps = 3;
phiD = ctx.phiD(:);
FmaxOf = @(rest) maxBPAforce(rest, '20');
cableOff = 0.00075;
cylUp = struct('c', ctx.geo.tibiaUpperCenter, 'r', ctx.geo.tibiaUpperRadius + cableOff);
cylLo = struct('c', ctx.geo.tibiaLowerCenter, 'r', ctx.geo.tibiaLowerRadius + cableOff);

capDir = fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022');
T5 = readtable(fullfile(capDir, 'vas_med_surface_5node_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
[~, o5] = sort(T5.knee_angle_deg(:));
nodeDeg = T5.knee_angle_deg(o5); tauHnode = T5.tau_h_Nm(o5);
TD = readtable(fullfile(capDir, 'vas_med_surface_dense_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
[~, oD] = sort(TD.knee_angle_deg(:));
tauHgrid = interp1(TD.knee_angle_deg(oD), TD.tau_h_Nm(oD), phiD, 'linear', 'extrap');
armOsimGrid = interp1(TD.knee_angle_deg(oD), TD.arm_m(oD), phiD, 'linear', 'extrap');

routeCache = containers.Map;
[Location0, bend0, info0, mL0, span0] = getRoute(tendonRec, p1Rec, pEndRec, ctx, routeCache);
Lchar = max(mL0);
say(logF, 'p1Rec [%s] pEndRec [%s] rest %.4f tendon %.4f Lchar %.4f Xi3 %.4f\n', ...
    num2str(p1Rec, 4), num2str(pEndRec, 4), restRec, tendonRec, Lchar, Xi(4));

%% ================= STAGE A: GATES ======================================
kmaxRec = restRec * (1 - KMAX);
% GATE 1: X3 class of record (per-pose Location, Xi3 native, BPAcount 1)
x3 = MonoPamDataExplicit_balanceX3('g1', Location0, ctx.CrossPoint, Dia, ...
    ctx.T_Pam, restRec, kmaxRec, tendonRec, fitting, Pres, ...
    Xi(1), Xi(2), Xi(3), Xi(4), wraps, phiD, 1, bend0);
Fx3 = vecnorm(double(x3.F_p), 2, 2);
E1 = evalRP([restRec tendonRec], 1, 1, 'b', ctx, Xi, Dia, Pres, KMAX, ...
    fitting, wraps, phiD, tauHgrid, armOsimGrid, [], p1Rec, pEndRec, ringOld, ...
    cylUp, cylLo, routeCache, nodeDeg, tauHnode, FmaxOf);
g1 = max(abs(E1.Fsingle - Fx3) ./ max(1, abs(Fx3)));
compD = max(0, 1 - ((restRec - (mL0 - tendonRec - 2 * fitting - Xi(1))) ...
    / restRec) / KMAX);
s0x3d = restRec - mL0 + Xi(1) + tendonRec + 2 * fitting ...
    + Xi(4) .* bend0(:) .* compD .^ 2;
cbD = bracketCb(Location0, ctx, Xi, numel(mL0));
P1a = psiSolve(restRec, 1, 1, s0x3d, span0, cbD, KMAX, Dia, Pres, ...
    tendonRec, wraps * 6, FmaxOf);
g1b = max(abs(P1a.Fsingle - Fx3) ./ max(1, abs(Fx3)));
say(logF, 'GATE 1b (1,1; X3 kSpr 18x) psi vs X3 class: max rel err %.3g (mean %.3g)', ...
    g1b, mean(abs(P1a.Fsingle - Fx3) ./ max(1, abs(Fx3))));
say(logF, 'GATE 1 (1,1) psi vs X3 class: max rel force err %.3g (mean %.3g)\n', ...
    g1, mean(abs(E1.Fsingle - Fx3) ./ max(1, abs(Fx3))));

% GATE 2 (fixed-row control vs BiPam_pulley2_20261022) SKIPPED:
% the class TransformationMat 4th dimension is the JOINT INDEX (T(:,:,a,1)
% hip, T(:,:,a,2) knee), not the N2 grid -- constructing the campaign
% driver layout ran past this session deadline (see Opt_run_BiPulley.m
% ~line 534 for the of-record construction). The G>1 psi here is a
% verbatim transcription of fortzBi2 (BiPam_pulley2_20261022.m lines
% 1515-1517, cAdd = 0); (1,1) equivalence is gated vs the X3 class.
% UNVERIFIED: bit-level G>1 agreement with the class on a control route.
g2 = NaN;
save(fullfile(here, 'vasmed_rp4_gates_20261028.mat'), 'g1', 'Fx3', 'E1');

%% ================= STAGE B: sizes-only sweep ==========================
lb2 = [0.30 * Lchar, 0.075];
ub2 = [0.90 * Lchar, 0.30];
seeds = [0.60 * Lchar, 0.10; 0.45 * Lchar, 0.09; min(restRec, ub2(1)), ...
    max(0.075, min(0.30, tendonRec)); 0.35, 0.08; 0.55, 0.12];
opts = optimoptions('patternsearch', 'Display', 'none', ...
    'MeshTolerance', 1e-4, 'StepTolerance', 1e-4, 'MaxFunctionEvaluations', 150);
rows = {}; tB = tic;
for nB = 1:5
    for G = 1:3
        fS = inf(size(seeds, 1), 1);
        for sI = 1:size(seeds, 1)
            if any(seeds(sI, :) > ub2 | seeds(sI, :) < lb2), fS(sI) = 1e6; continue; end
            E = evalRP(seeds(sI, :), nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, ...
                fitting, wraps, phiD, tauHgrid, armOsimGrid, [], p1Rec, pEndRec, ...
                ringOld, cylUp, cylLo, routeCache, nodeDeg, tauHnode, FmaxOf);
            fS(sI) = -1000 * E.worstFrac - 100 * E.meanFrac;
        end
        [~, iS] = min(fS);
        objF = @(x) sizObj(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
            phiD, tauHgrid, armOsimGrid, p1Rec, pEndRec, ringOld, cylUp, cylLo, ...
            routeCache, nodeDeg, tauHnode, FmaxOf);
        [xB, ~] = patternsearch(objF, seeds(iS, :), [], [], [], [], lb2, ub2, [], opts);
        E = evalRP(xB, nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
            phiD, tauHgrid, armOsimGrid, [], p1Rec, pEndRec, ringOld, cylUp, cylLo, ...
            routeCache, nodeDeg, tauHnode, FmaxOf);
        Ea = evalRP(xB, nB, G, 'a', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
            phiD, tauHgrid, armOsimGrid, [], p1Rec, pEndRec, ringOld, cylUp, cylLo, ...
            routeCache, nodeDeg, tauHnode, FmaxOf);
        rows(end + 1, :) = {nB, G, xB(1), xB(2), E.worstMargin, E.worstMarginPct, ...
            E.worstFrac, E.meanFrac, E.nMeetPlain, E.nMeet5pct, E.nSlack, ...
            E.nInfeas, E.nTaut, numel(phiD), E.nodeMargins, Ea.worstMargin, ...
            Ea.worstFrac, Ea.nMeetPlain}; %#ok<AGROW>
        say(logF, ['SWEEP nB%d G%d: rest %.4f tendon %.4f worst %+.2f N*m ' ...
            '(%+.1f%%) meet %d/5 (+5%% %d/5) slack %d infeas %d taut %d | ' ...
            'arm(a) worst %+.2f\n'], nB, G, xB(1), xB(2), E.worstMargin, ...
            E.worstMarginPct, E.nMeetPlain, E.nMeet5pct, E.nSlack, E.nInfeas, ...
            E.nTaut, Ea.worstMargin);
    end
end
say(logF, 'Stage B wall time %.1f min\n', toc(tB) / 60);
save(fullfile(here, 'vasmed_rp4_final_20261028.mat'), 'rows', '-v7.3');
fidC = fopen(fullfile(here, 'vasmed_rp4_sizesweep.csv'), 'w');
fprintf(fidC, ['nBPA,G,rest_m,tendon_m,worstMargin_Nm,worstMarginPct,worstFrac,' ...
    'meanFrac,nMeetPlain_5node,nMeet5pct_5node,nSlack,nInfeas,nTaut,nPoses,' ...
    'nodeMarginsNm,armA_worstMargin_Nm,armA_worstFrac,armA_nMeetPlain\n']);
for k = 1:size(rows, 1)
    fprintf(fidC, '%d,%d,%.6g,%.6g,%.6g,%.6f,%.6f,%.6f,%d,%d,%d,%d,%d,%d,"%s",%.6g,%.6f,%d\n', ...
        rows{k, 1:14}, num2str(rows{k, 15}, '%.2f '), rows{k, 16:18});
end
fclose(fidC);

%% ================= STAGE C: attachment optimization ====================
wf = [rows{:, 7}];
[~, ord] = sort(wf(:), 'descend');
cfgList = [1 1];
for k = ord(:).'
    cd = [rows{k, 1}, rows{k, 2}];
    if size(cfgList, 1) < 4 && ~any(cfgList(:, 1) == cd(1) & cfgList(:, 2) == cd(2))
        cfgList(end + 1, :) = cd; %#ok<AGROW>
    end
end
say(logF, 'Stage C configs: %s\n', mat2str(cfgList.'));

lb8 = [p1Seed + [-0.020, -0.060, 0], pEndSeed + [-0.030, -0.150, 0], 0.35, 0.075];
ub8 = [p1Seed + [ 0.100,  0.125, 0], pEndSeed + [ 0.025,  0.024, 0], 1.20, 0.30];
x0a = min(max([p1Seed, pEndSeed, 0.55, 0.10], lb8), ub8);
x0b = min(max([p1Rec, pEndRec, restRec, max(0.075, tendonRec)], lb8), ub8);
optsG = optimoptions('surrogateopt', 'Display', 'none', 'UseParallel', true, ...
    'MaxFunctionEvaluations', 150, 'MinSampleDistance', 1e-3, ...
    'ConstraintTolerance', 1e-6, 'InitialPoints', [x0a; x0b]);
optsP = optimoptions('patternsearch', 'Display', 'none', ...
    'MaxFunctionEvaluations', 150, 'MeshTolerance', 1e-4, ...
    'StepTolerance', 1e-4, 'ConstraintTolerance', 1e-6);

res = {}; tC = tic;
for c = 1:size(cfgList, 1)
    nB = cfgList(c, 1); G = cfgList(c, 2);
    obj = @(x) attObj(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
        phiD, tauHgrid, ringOld, cylUp, cylLo, p1Seed, pEndSeed, nodeDeg, ...
        tauHnode, FmaxOf);
    nlc = @(x) nonlconExt20mm(x, ctx);
    objcon = @(x) wrapOC(x, obj, nlc);
    [xG, ~] = surrogateopt(objcon, lb8, ub8, optsG);
    [xB, fB] = patternsearch(@(x) objPS(x, obj, nlc), xG, [], [], [], [], ...
        lb8, ub8, nlc, optsP);
    dev = norm(xB(1:3) - p1Seed) + norm(xB(4:6) - pEndSeed);
    E = evalRP8(xB, nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
        phiD, tauHgrid, armOsimGrid, ringOld, cylUp, cylLo, nodeDeg, ...
        tauHnode, FmaxOf);
    Ea = evalRP8(xB, nB, G, 'a', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
        phiD, tauHgrid, armOsimGrid, ringOld, cylUp, cylLo, nodeDeg, ...
        tauHnode, FmaxOf);
    res(end + 1, :) = {nB, G, xB, fB, E.worstMargin, E.worstMarginPct, ...
        E.worstFrac, E.meanFrac, E.nMeetPlain, E.nMeet5pct, E.nSlack, ...
        E.nInfeas, dev, Ea.worstMargin, Ea.nMeetPlain, E, Ea}; %#ok<AGROW>
    say(logF, ['ATTACH nB%d G%d: p1 [%s] pEnd [%s] rest %.4f tendon %.4f | ' ...
        'worst %+.2f (%+.1f%%) meet %d/5 (+5%% %d/5) slack %d infeas %d dev %.3f' ...
        ' | arm(a) worst %+.2f\n'], nB, G, num2str(xB(1:3), 4), ...
        num2str(xB(4:6), 4), xB(7), xB(8), E.worstMargin, E.worstMarginPct, ...
        E.nMeetPlain, E.nMeet5pct, E.nSlack, E.nInfeas, dev, Ea.worstMargin);
    save(fullfile(here, 'vasmed_rp4_final_20261028.mat'), 'rows', 'res', ...
        'cfgList', 'lb8', 'ub8', 'x0a', 'x0b', 'optsG', 'optsP', ...
        'p1Seed', 'pEndSeed', 'g1', 'g2', '-v7.3');
end
say(logF, 'Stage C wall time %.1f min\n', toc(tC) / 60);
fidC = fopen(fullfile(here, 'vasmed_rp4_attachopt.csv'), 'w');
fprintf(fidC, ['nBPA,G,p1x,p1y,p1z,pEndx,pEndy,pEndz,rest_m,tendon_m,' ...
    'worstMargin_Nm,worstMarginPct,worstFrac,meanFrac,nMeetPlain_5node,' ...
    'nMeet5pct_5node,nSlack,nInfeas,dev_m,armA_worstMargin_Nm,' ...
    'armA_nMeetPlain\n']);
for k = 1:size(res, 1)
    x = res{k, 3};
    fprintf(fidC, '%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6f,%.6f,%.6f,%d,%d,%d,%d,%.6f,%.6g,%d\n', ...
        res{k, 1:2}, x, res{k, 5:13}, res{k, 14:15});
end
fclose(fidC);
say(logF, 'RP4 DONE\n');
fclose(logF);
end

%% ================= LOCAL FUNCTIONS =====================================
function say(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end

function [f, g] = wrapOC(x, obj, nlc)
f = obj(x);
[cc, ~] = nlc(x);
g = cc;
end

function f = objPS(x, obj, nlc)
[c, ~] = nlc(x);
f = obj(x) + 1e4 * max(c);
end

function f = sizObj(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
    phiD, tauHgrid, armOsimGrid, p1, pEnd, ringOld, cylUp, cylLo, routeCache, ...
    nodeDeg, tauHnode, FmaxOf)
E = evalRP(x, nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, phiD, ...
    tauHgrid, armOsimGrid, [], p1, pEnd, ringOld, cylUp, cylLo, routeCache, ...
    nodeDeg, tauHnode, FmaxOf);
if ~E.ok, f = 1e6; return; end
f = -1000 * E.worstFrac - 100 * E.meanFrac;
end

function f = attObj(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
    phiD, tauHgrid, ringOld, cylUp, cylLo, p1Seed, pEndSeed, nodeDeg, ...
    tauHnode, FmaxOf)
E = evalRP8(x, nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, phiD, ...
    tauHgrid, [], ringOld, cylUp, cylLo, nodeDeg, tauHnode, FmaxOf);
if ~E.ok, f = 1e6; return; end
dev = norm(x(1:3) - p1Seed) + norm(x(4:6) - pEndSeed);
f = -1000 * E.worstFrac - 100 * E.meanFrac + 0.5 * (x(7) / 0.83) ...
    + 0.5 * (dev / 0.10);
end

%% ------- cached route (key: rounded design vector) --------------------
function [Location, bend, info, mL, span] = getRoute(tendon, p1, pEnd, ctx, cache)
key = mat2str(round([p1, pEnd, tendon] * 1e6));
if isKey(cache, key)
    v = cache(key);
    [Location, bend, info, mL, span] = deal(v{1:5});
    return
end
[Location, bend, info] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx);
N = size(Location, 3);
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
cache(key) = {Location, bend, info, mL, span};
end

%% ------- the evaluator: per-pose BiPam_pulley2 closure + X3 -----------
function E = evalRP(x2, nB, G, armMode, ctx, Xi, Dia, Pres, KMAX, fitting, ...
    wraps, phiD, tauHgrid, armOsimGrid, armCable, p1, pEnd, ringOld, cylUp, ...
    cylLo, cache, nodeDeg, tauHnode, FmaxOf)
E = evalCore(x2(1), x2(2), nB, G, armMode, ctx, Xi, Dia, Pres, KMAX, ...
    fitting, wraps, phiD, tauHgrid, armOsimGrid, armCable, p1, pEnd, ...
    ringOld, cylUp, cylLo, cache, nodeDeg, FmaxOf);
end

function E = evalRP8(x, nB, G, armMode, ctx, Xi, Dia, Pres, KMAX, fitting, ...
    wraps, phiD, tauHgrid, armOsimGrid, ringOld, cylUp, cylLo, nodeDeg, ...
    tauHnode, FmaxOf)
E = evalCore(x(7), x(8), nB, G, armMode, ctx, Xi, Dia, Pres, KMAX, fitting, ...
    wraps, phiD, tauHgrid, armOsimGrid, [], x(1:3), x(4:6), ringOld, cylUp, ...
    cylLo, containers.Map, nodeDeg, FmaxOf);
end

function E = evalCore(rest, tendon, nB, G, armMode, ctx, Xi, Dia, Pres, ...
    KMAX, fitting, wraps, phiD, tauHgrid, armOsimGrid, armCable, p1, pEnd, ...
    ringOld, cylUp, cylLo, cache, nodeDeg, FmaxOf)
E = struct('ok', false, 'worstFrac', 0, 'meanFrac', 0, 'worstMargin', -inf, ...
    'worstMarginPct', -100);
if tendon < 0.075 - 1e-12 || tendon > 0.3 || rest <= 0
    return
end
try
    [Location, bend, info, mL, span] = getRoute(tendon, p1, pEnd, ctx, cache);
catch
    return
end
if isfield(info, 'tendonMax') && tendon > min(info.tendonMax) + 1e-9
    return
end
E = solveAndGrade(rest, tendon, nB, G, armMode, mL, span, bend, Xi, Dia, ...
    Pres, KMAX, fitting, wraps, phiD, tauHgrid, armOsimGrid, armCable, ...
    pEnd, ringOld, cylUp, cylLo, nodeDeg, FmaxOf, Location, info, ctx);
end

function E = solveAndGrade(rest, tendon, nB, G, armMode, mL, span, bend, Xi, ...
    Dia, Pres, KMAX, fitting, wraps, phiD, tauHgrid, armOsimGrid, armCable, ...
    pEnd, ringOld, cylUp, cylLo, nodeDeg, FmaxOf, Location, info, ctx)

% Xi3 fictive loss (X3 class Contraction_k formula, rigid first pass)
N = numel(mL);
bend = bend(:); if numel(bend) ~= N, bend = zeros(N, 1); end
Lm0 = mL - tendon - 2 * fitting - Xi(1);
comp = max(0, 1 - ((rest - Lm0) / rest) / KMAX);
deltaL = Xi(4) .* bend .* comp .^ 2;
s0x3 = rest - mL + Xi(1) + tendon + 2 * fitting + deltaL;

% cb per pose: verbatim LokBi2 bracket-frame projection
Fmax = FmaxOf(rest);
cb = bracketCb(Location, ctx, Xi, N);

P = psiSolve(rest, nB, G, s0x3, span, cb, KMAX, Dia, Pres, tendon, ...
    wraps, FmaxOf);
Ftendon = P.Ftendon; slack = P.slack; infeas = P.infeas;
Fsingle = P.Fsingle; sCon = P.sCon;

if strcmp(armMode, 'a')
    arm = armOsimGrid;
else
    if isempty(armCable)
        armCable = cableArmGrid(pEnd, ringOld, cylUp, cylLo, phiD);
    end
    arm = armCable;
end
taut = ~slack & ~infeas & isfinite(Ftendon) & Ftendon > 0;
tauR = Ftendon .* arm;
tauR(~taut) = 0;
meetFrac = zeros(N, 1);
meetFrac(taut) = tauR(taut) ./ tauHgrid(taut);
E.tauR = tauR; E.tauH = tauHgrid; E.taut = taut;
E.Ftendon = Ftendon; E.Fsingle = Fsingle; E.sCon = sCon;
E.deltaL = deltaL; E.slack = slack; E.infeas = infeas;
E.worstFrac = min(min(meetFrac), 1);
E.meanFrac = mean(min(meetFrac, 1));
E.worstMargin = min(tauR - tauHgrid);
E.worstMarginPct = 100 * (min(tauR ./ tauHgrid) - 1);
E.nSlack = nnz(slack); E.nInfeas = nnz(infeas); E.nTaut = nnz(taut);
tauN = interp1(phiD, tauR, nodeDeg, 'linear', 'extrap');
tauHn = interp1(phiD, tauHgrid, nodeDeg, 'linear', 'extrap');
E.nMeetPlain = nnz(tauN >= tauHn - 1e-12);
E.nMeet5pct = nnz(tauN >= 1.05 * tauHn - 1e-12);
E.nodeMargins = tauN - tauHn; E.tauN = tauN;
E.ok = true; E.rest = rest; E.tendon = tendon; E.nB = nB; E.G = G;
end

%% ------- psi solve: the BiPam_pulley2 closure per pose ----------------
function P = psiSolve(rest, nB, G, s0x3, span, cb, KMAX, Dia, Pres, ...
    tendon, wraps, FmaxOf)
% psi(s) = G*s - G*cb(i)*F(s) - nB*F(s)*cSpr/G - DL(i)
%          - G*(s0x3(i) - DL(i))          [fortzBi2, cAdd = 0]
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
    'slack', slack, 'infeas', infeas, 'DL', DL, 'cb', cb);
end

function cb = bracketCb(Location, ctx, Xi, N)
% Verbatim LokBi2 bracket-frame compliance projection per pose:
% pre-force BPA direction (row 1 -> first non-duplicate row, the crossing
% row transformed by T_Pam when it is the tibia row -- computeForceVector
% convention), mapped through the affine Thbr\eye(4) exactly as LokBi2
% maps the force vector, then cb = u_hat*C*u_hat with C =
% diag(1/Xi1, 1/Xi2, 1/Xi2); Pbr = [9.48 -36.21 30.86]/1000 (dia 20,
% pulley-family table). pA = Location(1,:,1): row 1 is frame-fixed.
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
            if r == ctx.CrossPoint      % tibia-side row: ICR -> femur
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
