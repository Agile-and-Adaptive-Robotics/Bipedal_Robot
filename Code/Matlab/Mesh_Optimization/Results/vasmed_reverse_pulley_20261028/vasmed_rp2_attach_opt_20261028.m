function vasmed_rp2_attach_opt_20261028()
% vasmed_rp2_attach_opt_20261028 -- STUDY (2026-10-28, easteregg2)
%
% Part 2 of the reverse-pulley vastus medialis example: ATTACHMENT-PLACEMENT
% OPTIMIZATION. The 10-27 study1 (regraded in rp1) re-solved only
% [rest tendon] on the FIXED route of record -- a length-sizing sweep,
% not placement optimization. This script optimizes the permitted
% attachments on the canonical route machinery:
%
%   x = [p1(1:3), pEnd(1:3), rest, tendon]     (8 dimensions)
%
% p1 = femur-frame BPA ORIGIN (route row 1); pEnd = t1-frame tibial
% INSERTION (route row p9) -- both consumed by the route builder of
% record buildDistalRingLocation20mm. This is the TIBIA-ROUTE solver:
% its insertion is a free tibial attachment, so both endpoints are
% optimizable (Ben, 2026-10-28). The OTHER solver family (Routing27 /
% OpenSim robot-body paths whose insertion is the patella, given as a
% femur-frame function of knee angle) cannot optimize its insertion and
% is not used here.
%
% Bounds (documented sources, buildKneeExtContext20mm.m):
%   p1    : pre-pin p1 search region p1_0 + [-0.020 -0.060 0] ..
%           [+0.100 +0.125 0], p1_0 = routeSeed row 1. The cam campaign
%           later pinned p1 = [0.040 0.137 0]; Ben's 2026-10-28 ruling
%           that the resulting 0.8277 m rest is too long re-opens the
%           documented region so the origin may move distal.
%   pEnd  : active bounds of record pEnd_0 + [-0.030 -0.150 0] ..
%           [+0.025 +0.024 0].
%   rest  : [0.35, 1.2] m with penalty 0.5*rest/0.83 (shorter wins ties).
%   tendon: [0.075, 0.3] m. 0.075 m = pulley packaging floor (2026-10-28
%           directive); it is the ROUTE-MODEL tendon parameter only --
%           the physical line looped through a G-part tackle and
%           terminated at the insertion is longer by the tackle's
%           internal line, which the route model does not carry.
%
% Deviation penalty (Ben, 2026-10-28): +0.5*dev/0.10,
% dev = ||p1-p1_seed|| + ||pEnd-pEnd_seed|| against the anatomical seed
% rows (routeSeed 1 and 9). Weight small: a 0.01 worst-domain meet gain
% (10 units) dominates any plausible deviation/rest change.
%
% Constraints: nonlconExt20mm(x, ctxCamBest) -- the canonical constraint
% vector (geometry-dependent tendon max, minimum-possible tendon, pEnd
% envelope, p8 x-limit, zero-strain rest feasibility, Contraction/KMAX
% <= 1); it evaluates MonoPamDataExplicit_balanceX3 with the adopted Xi3
% natively through predictKneeExt20mm.
%
% Objective (arm method b, routed cable arm about the moving ICR):
%   f = -1000*worstFrac - 100*meanFrac + 0.5*(rest/0.83) + 0.5*(dev/0.10)
% graded over ALL 100 poses of -120..+10 deg, X3-corrected force,
% slack/infeasible/zero-force poses = full shortfalls.
% Solver: surrogateopt (250 evals, parallel, MinSampleDistance 1e-3,
% ConstraintTolerance 1e-6) + patternsearch polish (250 evals,
% MeshTolerance/StepTolerance 1e-4) -- the route-redesign solver pair
% of record.
%
% Configs: top 3 (nB,G) cells by rp1 arm(b) worstFrac + direct-drive
% (1,1), from vasmed_rp1_regrade_20261028.mat.
%
% Outputs: vasmed_rp2_attachopt_20261028.mat / ..._summary.csv / log.
%
% Invocation: matlab -batch "cd(<folder>); vasmed_rp2_attach_opt_20261028"

R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
here = fileparts(mfilename('fullpath'));
logF = fopen(fullfile(here, 'vasmed_rp2_attachopt_20261028.txt'), 'w');

M1 = load(fullfile(here, 'vasmed_rp1_regrade_20261028.mat'));
rows = M1.rows;
worstB = [rows{:, 7}];                    % arm(b) worstFrac per config
nBv = [rows{:, 1}]; Gv = [rows{:, 2}];
[~, ord] = sort(worstB(:), 'descend');
cfgList = [1 1];                          % always the direct-drive cell
for k = ord(:).'
    cand = [nBv(k), Gv(k)];
    if size(cfgList, 1) < 4 && ~any(cfgList(:, 1) == cand(1) & cfgList(:, 2) == cand(2))
        cfgList(end + 1, :) = cand; %#ok<AGROW>
    end
end
say(logF, 'Configs to optimize (nB G): %s\n', mat2str(cfgList.'));

% ---- fixed data ----
C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'ctxCamBest', 'XiUsed');
ctx = C.ctxCamBest;
Xi = C.XiUsed(:).';
xRec = C.xFullBest(:).';
p1Rec = xRec(1:3); pEndRec = xRec(4:6);
ringOld = ctx.routeSeed(6, :);
p1Seed = ctx.routeSeed(1, :);
pEndSeed = ctx.routeSeed(9, :);
Dia = 20; Pres = 620; KMAX = 0.255; fitting = 0.021; wraps = 3;
phiD = ctx.phiD(:);
tauHgrid = M1.tauHgrid;
armOsimGrid = interp1(M1.denseDeg, M1.armOsimDense, phiD, 'linear', 'extrap');
nodeDeg = M1.nodeDeg; tauHnode = M1.tauHnode;
FmaxOf = @(rest) maxBPAforce(rest, '20');
cableOff = 0.00075;
cylUp = struct('c', ctx.geo.tibiaUpperCenter, 'r', ctx.geo.tibiaUpperRadius + cableOff);
cylLo = struct('c', ctx.geo.tibiaLowerCenter, 'r', ctx.geo.tibiaLowerRadius + cableOff);

% ---- bounds ----
lb = [p1Seed + [-0.020, -0.060, 0], pEndSeed + [-0.030, -0.150, 0], 0.35, 0.075];
ub = [p1Seed + [ 0.100,  0.125, 0], pEndSeed + [ 0.025,  0.024, 0], 1.20, 0.30];
x0 = [p1Seed, pEndSeed, max(0.35, min(1.2, M1.restRec)), ...
    max(0.075, min(0.30, M1.tendonRec))];
x0 = min(max(x0, lb), ub);
xFixed = min(max([p1Rec, pEndRec, M1.restRec, M1.tendonRec], lb), ub);
say(logF, ['lb %s\nub %s\nseeds: anatomical %s ; design-of-record(clipped) %s\n'], ...
    mat2str(lb, 4), mat2str(ub, 4), mat2str(x0, 4), mat2str(xFixed, 4));

% ---- solver options of record ----
optsG = optimoptions('surrogateopt', 'Display', 'none', 'UseParallel', true, ...
    'MaxFunctionEvaluations', 250, 'MinSampleDistance', 1e-3, ...
    'ConstraintTolerance', 1e-6, 'InitialPoints', [x0; xFixed]);
optsP = optimoptions('patternsearch', 'Display', 'none', ...
    'MaxFunctionEvaluations', 250, 'MeshTolerance', 1e-4, ...
    'StepTolerance', 1e-4, 'ConstraintTolerance', 1e-6);

res = {}; tAll = tic;
for c = 1:size(cfgList, 1)
    nB = cfgList(c, 1); G = cfgList(c, 2);
    obj = @(x) attachObj(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
        phiD, tauHgrid, ringOld, cylUp, cylLo, p1Seed, pEndSeed, nodeDeg, ...
        tauHnode, FmaxOf);
    nlc = @(x) nonlconExt20mm(x, ctx);
    objcon = @(x) wrapObjCon(x, obj, nlc);
    [xG, ~] = surrogateopt(objcon, lb, ub, optsG);
    [xB, fB] = patternsearch(@(x) objPS(x, obj, nlc), xG, [], [], [], [], ...
        lb, ub, nlc, optsP);
    dev = norm(xB(1:3) - p1Seed) + norm(xB(4:6) - pEndSeed);
    E = evalAttach(xB, nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
        phiD, tauHgrid, armOsimGrid, ringOld, cylUp, cylLo, nodeDeg, ...
        tauHnode, FmaxOf);
    Ea = evalAttach(xB, nB, G, 'a', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
        phiD, tauHgrid, armOsimGrid, ringOld, cylUp, cylLo, nodeDeg, ...
        tauHnode, FmaxOf);
    res(end + 1, :) = {nB, G, xB, fB, E.worstMargin, E.worstMarginPct, ...
        E.worstFrac, E.meanFrac, E.nMeetPlain, E.nMeet5pct, E.nSlack, ...
        E.nInfeas, dev, Ea.worstMargin, Ea.nMeetPlain, E, Ea}; %#ok<AGROW>
    say(logF, ['ATTACH nB%d G%d: p1 [%s] pEnd [%s] rest %.4f tendon %.4f | ' ...
        'worst %+.2f N*m (%+.1f%%) meet %d/%d (+5%% %d/%d) slack %d infeas %d ' ...
        'dev %.3f m | arm(a) worst %+.2f\n'], nB, G, num2str(xB(1:3), 4), ...
        num2str(xB(4:6), 4), xB(7), xB(8), E.worstMargin, E.worstMarginPct, ...
        E.nMeetPlain, numel(nodeDeg), E.nMeet5pct, E.nSlack, E.nInfeas, dev, ...
        Ea.worstMargin);
    save(fullfile(here, 'vasmed_rp2_attachopt_20261028.mat'), 'res', 'cfgList', ...
        'lb', 'ub', 'x0', 'xFixed', 'optsG', 'optsP', 'p1Seed', 'pEndSeed');
end
say(logF, 'Attach-opt wall time: %.1f min\n', toc(tAll) / 60);

fidC = fopen(fullfile(here, 'vasmed_rp2_attachopt_summary_20261028.csv'), 'w');
fprintf(fidC, ['nBPA,G,p1x,p1y,p1z,pEndx,pEndy,pEndz,rest_m,tendon_m,' ...
    'worstMargin_Nm,worstMarginPct,worstFrac,meanFrac,nMeetPlain_5node,' ...
    'nMeet5pct_5node,nSlack,nInfeas,dev_m,armA_worstMargin_Nm,' ...
    'armA_nMeetPlain\n']);
for k = 1:size(res, 1)
    x = res{k, 3};
    fprintf(fidC, '%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6f,%.6f,%.6f,%d,%d,%d,%d,%.6f,%.6g,%d\n', ...
        res{k, 1}, res{k, 2}, x, res{k, 5:14}, res{k, 15:16});
end
fclose(fidC);
say(logF, 'RP2 DONE\n');
fclose(logF);
end

%% ================= LOCAL FUNCTIONS =====================================
function say(fid, fmt, varargin)
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end

function [f, g] = wrapObjCon(x, obj, nlc)
% surrogateopt objective-with-constraints form: g <= 0 feasible.
f = obj(x);
[cc, ~] = nlc(x);
g = cc;
end

function f = objPS(x, obj, nlc)
% patternsearch objective: hard-penalize constraint violation.
[c, ~] = nlc(x);
f = obj(x) + 1e4 * max(c);
end

function f = attachObj(x, nB, G, ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
    phiD, tauHgrid, ringOld, cylUp, cylLo, p1Seed, pEndSeed, nodeDeg, ...
    tauHnode, FmaxOf)
E = evalAttach(x, nB, G, 'b', ctx, Xi, Dia, Pres, KMAX, fitting, wraps, ...
    phiD, tauHgrid, [], ringOld, cylUp, cylLo, nodeDeg, tauHnode, FmaxOf);
if ~E.ok
    f = 1e6; return
end
dev = norm(x(1:3) - p1Seed) + norm(x(4:6) - pEndSeed);
f = -1000 * E.worstFrac - 100 * E.meanFrac + 0.5 * (x(7) / 0.83) ...
    + 0.5 * (dev / 0.10);
end

%% ------- corrected equilibrium evaluator, ATTACHMENT-full design ------
function E = evalAttach(x, nB, G, armMode, ctx, Xi, Dia, Pres, KMAX, ...
    fitting, wraps, phiD, tauHgrid, armOsimGrid, ringOld, cylUp, cylLo, ...
    nodeDeg, tauHnode, FmaxOf)
p1 = x(1:3); pEnd = x(4:6);
rest = x(7); tendon = x(8);
E = struct('ok', false, 'worstFrac', 0, 'meanFrac', 0, 'worstMargin', -inf, ...
    'worstMarginPct', -100);
if tendon < 0.075 - 1e-12 || tendon > 0.3 || rest <= 0
    return
end
try
    [Location, bendM, info] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx);
catch
    return
end
if isfield(info, 'tendonMax') && tendon > min(info.tendonMax) + 1e-9
    return
end
mL = zeros(size(phiD, 1), 1);
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
obj = MonoPam_pulley('rp2', Location, ctx.CrossPoint, Dia, ctx.T_Pam, ...
    rest, rest * (1 - KMAX), tendon, fitting, Pres, Xi(1), Xi(2), Xi(3), ...
    wraps, pulCfg);
sCon = double(obj.sContraction(:));
slack = logical(obj.PulleySlack(:));
infeas = logical(obj.PulleyInfeasible(:));
bend = bendM(:);
if numel(bend) ~= N, bend = zeros(N, 1); end
Lm0 = mL - tendon - 2 * fitting - Xi(1);
strain0 = (rest - Lm0) / rest;
comp = max(0, 1 - (strain0 / KMAX));
deltaL = Xi(4) .* bend .* comp .^ 2;
Fsingle = max(festo4(Dia, (sCon + deltaL) / rest / KMAX, Pres) .* FmaxOf(rest), 0);
Ftendon = nB * Fsingle / G;
Ftendon(slack) = 0;
Ftendon(infeas) = NaN;

if strcmp(armMode, 'a')
    arm = armOsimGrid;
else
    arm = cableArmGrid(pEnd, ringOld, cylUp, cylLo, phiD);
end
taut = ~slack & ~infeas & isfinite(Ftendon) & Ftendon > 0;
tauR = Ftendon .* arm;
tauR(~taut) = 0;
meetFrac = zeros(N, 1);
meetFrac(taut) = tauR(taut) ./ tauHgrid(taut);
E.tauR = tauR; E.tauH = tauHgrid; E.taut = taut;
E.Ftendon = Ftendon; E.arm = arm; E.deltaL = deltaL;
E.worstFrac = min(min(meetFrac), 1);
E.meanFrac = mean(min(meetFrac, 1));
E.worstMargin = min(tauR - tauHgrid);
E.worstMarginPct = 100 * (min(tauR ./ tauHgrid) - 1);
E.nSlack = nnz(slack); E.nInfeas = nnz(infeas); E.nTaut = nnz(taut);
E.slack = slack; E.infeas = infeas;
tauN = interp1(phiD, tauR, nodeDeg, 'linear', 'extrap');
tauHn = interp1(phiD, tauHgrid, nodeDeg, 'linear', 'extrap');
meetP = tauN >= tauHn - 1e-12;
meet5 = tauN >= 1.05 * tauHn - 1e-12;
E.nMeetPlain = nnz(meetP); E.nMeet5pct = nnz(meet5);
E.nodeMargins = tauN - tauHn; E.tauN = tauN;
E.ok = true; E.x = x;
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
