function [f, rep] = fr_objective(C, cs, kin, x, xiBranch, opts)
%FR_OBJECTIVE Spec-driven objective + hard constraints for one candidate.
%
% The objective TERMS AND WEIGHTS are read FRESH from objective_spec.json
% at every solve (master item 8: never hardcoded, iterable). Every
% evaluation PRINTS the terms separately with their weights and
% normalizations (rep.print / console line when ~opts.quiet).
%
% PRIMARY objectives (master item 8): torque capacity (meet-or-beat the
% human target torque on the primary dof) and usable motion range. A
% small, separately weighted torque-shape penalty is defined apart from
% the capacity terms. BPA-BPA overlap is a SOFT penalty from the
% pose-dependent inflated envelopes (physically simultaneous poses only).
%
% FORBIDDEN and structurally absent: any off-axis torque penalty (master
% item 8 / handoff item 9). Secondary joint torques are RECORDED in
% rep.tauSecondary and never enter f - the objective-spec file lists this
% as forbidden, and no term of that class exists in the registry.
%
% HARD constraints (violation form, c <= 0 feasible):
%   1 z=0 pelvis-frame centerline (handoff item 8, every attachment and
%     mount anchor at every evaluated pose, sideSign*z >= 0)
%   2 foot occupancy ban (sleeve fitting tip stays proximal of the ankle
%     pivot plane; BPA/tackle owning bodies structurally off the foot)
%   3 required tendon spans (fr_spans, computed from geometry + hardware)
%   4 transmission feasibility (no infeasible/slack cells; strain floor;
%     gain <= line parts structurally)
%   5 rated hardware loads (tackle vector reaction and the 180-deg
%     reversing-sheave equivalent 2T vs the assigned family's catalog
%     capacity)
%   6 bone exclusion (point-cloud distance >= inflated radius + rigid
%     clearance; buildGeoExclusion/nonlconExclusion lineage)
%   7 retained deeper-layer envelopes (Ben's layering ruling; HARD) with
%     the documented CAD-registration hook
%
% f = weighted objective + PENALTY * total constraint violation.

if nargin < 6
    opts = struct();
end
opts = fr_defaultOpts(opts);

tEval = tic;
% constraint accumulator: generously preallocated, trimmed at the end
% (max realistic count: 1 centerline + 1 foot + nPasses spans + 3 tx +
%  2 tackle + 2 obstacles << 32)
NC_MAX = 32;
cvec = nan(NC_MAX, 1);
clabels = cell(NC_MAX, 1);
ci = 0;
% ---- 0. objective spec (fresh read: iterable, never hardcoded) ----------
spec = jsondecode(fileread(fullfile(C.camp, 'objective_spec.json')));
termsList = fr_termList(spec.terms);   % jsondecode may give cell or struct
W = containers.Map();
for i = 1:numel(termsList)
    t = termsList{i};
    W(t.name) = t.weight;
end

% ---- 1. design + score-pass mechanics -----------------------------------
design = fr_unpackX(cs, x);
scorePass = 1;
if isfield(opts, 'passIdx') && ~isempty(opts.passIdx)
    scorePass = opts.passIdx;
elseif isfield(cs, 'scorePass') && ~isempty(cs.scorePass)
    scorePass = cs.scorePass;
end
qg = struct('passIdx', scorePass, 'qList', ...
    {cs.passes{scorePass}.qStructs});
[mech, meta] = fr_mechBuild(C, cs, design, kin, xiBranch, qg);
nG = numel(qg.qList);

% ---- 2. torques by virtual work on the solved state ----------------------
% tau_J(i) = segTension(i,J) * dDeltaL_J/dtheta (DeltaL = span(theta0) -
% span(theta), so dDeltaL/dtheta = -dSpan/dtheta and tau = -T*dSpan/dtheta,
% the mech4-verified chord virtual-work convention). The class's own
% Torque_ins is carried as a cross-check (reported).
kPrim = cs.primaryCrossing;
dTheta = cs.passes{scorePass}.th1(2) - cs.passes{scorePass}.th1(1);
% torque by CHORD virtual work: tau = -T * dChord/dtheta, the exact
% virtual-work partner per the mech4 verification (the polyline deltaL
% derivative is NOT the exact partner - documented mech4 'e-straight'
% expected mismatch). Chord = deformed exit->insertion distance in the
% crossing's child frame, from the class's pointInFrameDeformed.
chord = fr_crossingChord(mech, kPrim);
tauPrim = -mech.segTension(:, kPrim) .* fr_dDeltaL(chord, dTheta);
tauPrim = edgeFix(tauPrim);
tauClass = zeros(nG, 1);
for i = 1:nG
    % per-pose insertion torque about the primary pivot, from the class's
    % Torque_ins at grid cell (i, primary crossing)
    tq = mech.cell.Torque_ins(i, :, 1, kPrim);
    tauClass(i) = dot(reshape(tq, 1, 3), reshape(cs.primaryAxis, 1, 3)) ...
        / max(norm(cs.primaryAxis), eps);
end
tauSecondary = [];
if cs.secondaryCrossing > 0
    k2 = cs.secondaryCrossing;
    chord2 = fr_crossingChord(mech, k2);
    tauSecondary = edgeFix(-mech.segTension(:, k2) .* ...
        fr_dDeltaL(chord2, dTheta));
end

% ---- 3. objective terms (spec-driven) ------------------------------------
tauT = cs.target.tau(:);
if numel(tauT) ~= nG
    tauT = interpTarget(tauT, nG);
end
normCap = sum(abs(tauT));
deficitVec = max(0, abs(tauT) - abs(tauPrim));
termDeficit = sum(deficitVec) / max(normCap, eps);
termRom = mean(deficitVec > 1e-6);
termShape = mean((tauPrim - tauT).^2) / max(mean(tauT.^2), eps);
[termOverlap, overlapMaxMm] = fr_softOverlap(C, cs, design, kin, mech, meta);

terms = struct();
terms.torque_capacity_deficit = termDeficit;
terms.usable_rom_shortfall = termRom;
terms.torque_shape = termShape;
terms.bpa_bpa_overlap_soft = termOverlap;

fObj = 0;
names = fieldnames(terms);
for i = 1:numel(names)
    fObj = fObj + W(names{i}) * terms.(names{i});
end

% ---- 4. hard constraints -------------------------------------------------

% 1. centerline z=0 (every attachment + anchor, every pass pose, per-pass
% used-row indexing so masks of different passes each index correctly)
zmin = inf;
for p = 1:numel(cs.passes)
    qL = cs.passes{p}.qStructs;
    [locP, ~] = fr_locForPass(C, cs, design, kin, p);
    for i = 1:numel(qL)
        for u = 1:size(locP, 1)
            zw = fr_rowZPelvis(cs, kin, qL{i}, p, u, locP(u, :));
            zmin = min(zmin, zw);
        end
    end
end
ci = ci + 1; cvec(ci) = -zmin; clabels{ci} = 'centerline_z0'; 

% 2. foot occupancy (fitting-tip rule; owning bodies structural)
if cs.footTipRule
    tipMin = fr_fittingTipMargin(C, cs, kin, meta);
    ci = ci + 1; cvec(ci) = -tipMin; clabels{ci} = 'foot_occupancy_ban';
end

% 3. required tendon spans (master item 6): the OPTIMIZED routing must
% never need more tendon than the PROCURED hardware quantity
% cs.tendonHardware_m = max seed chord + ferrule/stop allowance (the
% contract placeholder resolved at case setup from the seed geometry)
for p = 1:numel(cs.passes)
    if ~strcmp(cs.spanKind, 'none')
        [ptsP, bodP] = fr_locForPass(C, cs, design, kin, p);
        sp = fr_spans(C, kin, ptsP, bodP, cs.passes{p}.qStructs, cs.spanKind);
        ci = ci + 1;
        cvec(ci) = (sp.spanMax_m + sp.allowance_m - ...
            cs.tendonHardware_m) / max(cs.tendonHardware_m, eps);
        clabels{ci} = sprintf('span_%s_pass%d', cs.spanKind, p);
        if p == 1
            rep.spanRequired = sp;
            rep.tendonHardware_m = cs.tendonHardware_m;
            rep.spanChordMax_m = sp.spanMax_m;
        end
    end
end

% 4. transmission feasibility over ALL passes
infeasFrac = 0; slackFrac = 0; elongViol = 0;
for p = 1:numel(cs.passes)
    qgp = struct('passIdx', p, 'qList', {cs.passes{p}.qStructs});
    [mp, ~] = fr_mechBuild(C, cs, design, kin, xiBranch, qgp);
    infeasFrac = max(infeasFrac, mean(mp.PulleyInfeasible));
    slackFrac = max(slackFrac, mean(mp.PulleySlack));
    elongViol = max(elongViol, max(0, -min(mp.strain_p) - strainFloor(cs)));
end
ci = ci + 1; cvec(ci) = infeasFrac; clabels{ci} = 'tx_infeasible_cells';
ci = ci + 1; cvec(ci) = slackFrac; clabels{ci} = 'tx_slack_cells';
ci = ci + 1; cvec(ci) = elongViol; clabels{ci} = 'strain_floor';

% 5. rated hardware loads
if cs.hasTackle
    fam = C.families(cs.txConfig.tackleFamily);
    cap = fam.capacity_N;
    loadMax = max(mech.ReactionFmag);
    revMax = 2 * max(mech.Ftendon);   % 180-deg reversing sheave worst case
    ci = ci + 1; cvec(ci) = (loadMax - cap) / cap;
    clabels{ci} = sprintf('tackle_load_%s', cs.txConfig.tackleFamily);
    ci = ci + 1; cvec(ci) = (revMax - cap) / cap;
    clabels{ci} = 'reversing_sheave_2T';
    rep.tackleLoadMax_N = loadMax;
    rep.reversingSheaveMax_N = revMax;
    rep.tackleCapacity_N = cap;
end

% 6 + 7. bone + retained envelopes (all passes' poses, physically
% simultaneous per pass)
obsBoneMin = inf; obsEnvMin = inf;
for p = 1:numel(cs.passes)
    [ptsP, bodP] = fr_locForPass(C, cs, design, kin, p);
    csP = cs; csP.cross = cs.crossFor{p};
    ob = fr_obstacles(C, csP, ptsP, bodP, kin, ...
        cs.passes{p}.qStructs, opts.envelopes);
    obsBoneMin = min(obsBoneMin, ob.boneMargin);
    obsEnvMin = min(obsEnvMin, ob.envMargin);
    if p == 1
        rep.obstacles = ob;
    end
end
ci = ci + 1; cvec(ci) = -obsBoneMin; clabels{ci} = 'bone_exclusion';
ci = ci + 1; cvec(ci) = -obsEnvMin; clabels{ci} = 'retained_envelopes';

cvec = cvec(1:ci); clabels = clabels(1:ci);
viol = max(0, cvec);
f = fObj + opts.penaltyScale * sum(viol);

% ---- 5. report ------------------------------------------------------------
rep = fr_reportShell(rep);
rep.terms = terms;
rep.weights = struct();
rep.normalization = struct();
for i = 1:numel(spec.terms)
    rep.weights.(termsList{i}.name) = termsList{i}.weight;
    rep.normalization.(termsList{i}.name) = termsList{i}.normalization;
end
rep.fObj = fObj;
rep.f = f;
rep.violation = sum(viol);
rep.c = cvec;
rep.clabels = clabels;
rep.tauPrimary = tauPrim;
rep.tauClassCrossCheck = tauClass;
rep.tauCrossCheckMaxDiff = max(abs(tauPrim - tauClass));
rep.tauTarget = tauT;
rep.tauSecondary = tauSecondary;
rep.meetsAtAllPoses = all(deficitVec <= 1e-6);
rep.fractionPosesMet = 1 - termRom;
rep.overlapMaxMm = overlapMaxMm;
rep.design = design;
rep.loc = meta.loc;
rep.wallSeconds = toc(tEval);
rep.specStatus = spec.status;
rep.midPoseWorld = fr_rowWorld(cs, kin, ...
    cs.passes{scorePass}.qStructs{ceil(nG / 2)}, scorePass, meta.loc);
rep.print = sprintf( ...
    ['EVAL %-16s t=%.3fs f=%+.4f (obj %+.4f viol %.2e) ' ...
    'deficit=%.3f(w%g,n=Sum|tauT|) rom=%.2f(w%g,n=fraction) ' ...
    'shape=%.3f(w%g,n=mean/mean) overlap=%.3f(w%g,n=Drest; max %.2f mm) ' ...
    'cmax=%.2e'], ...
    cs.name, rep.wallSeconds, f, fObj, sum(viol), ...
    termDeficit, W('torque_capacity_deficit'), ...
    termRom, W('usable_rom_shortfall'), ...
    termShape, W('torque_shape'), ...
    termOverlap, W('bpa_bpa_overlap_soft'), overlapMaxMm, max(viol));
if ~opts.quiet
    fprintf('%s\n', rep.print);
end
end

% ==========================================================================
function opts = fr_defaultOpts(opts)
if ~isfield(opts, 'envelopes'); opts.envelopes = {}; end
if ~isfield(opts, 'quiet'); opts.quiet = false; end
if ~isfield(opts, 'penaltyScale'); opts.penaltyScale = 1e3; end
end

function tl = fr_termList(terms)
% jsondecode returns a JSON object array as a cell of structs when the
% objects are heterogeneous (and as a struct array otherwise); normalize
tl = cell(1, numel(terms));
for i = 1:numel(terms)
    if iscell(terms)
        tl{i} = terms{i};
    else
        tl{i} = terms(i);
    end
end
end


function c = fr_crossingChord(mech, k)
% deformed exit->insertion chord per pose in crossing k's child frame
m = mech.cell;
exitRow = m.PulleyExitIndex(k);
insRow = size(m.Location, 1);
nG = mech.nG;
c = zeros(nG, 1);
for i = 1:nG
    pE = m.pointInFrameDeformed(exitRow, k, i, 1);
    pI = m.pointInFrameDeformed(insRow, k, i, 1);
    c(i) = norm(pI - pE);
end
end
function d = fr_dDeltaL(v, dTheta)
% d(v)/d(theta): v = deltaL (span change from the theta0 reference), so
% tau = -T*dSpan/dtheta = +T*dDeltaL/dtheta (the mech4-verified chord
% virtual-work convention). Central differences over the uniform grid.
n = numel(v);
d = zeros(n, 1);
d(2:n-1) = (v(3:n) - v(1:n-2)) / (2 * dTheta);
d(1) = (v(2) - v(1)) / dTheta;
d(n) = (v(n) - v(n-1)) / dTheta;
end

function v = edgeFix(v)
v = v(:);
end

function floorS = strainFloor(cs)
if strcmp(cs.diaStr, '10')
    floorS = 0.02;
else
    floorS = 0.03;
end
end

function [depth, depthMaxMm] = fr_softOverlap(C, cs, ~, ~, mech, meta)
% soft BPA-BPA overlap between the group's OWN parallel sleeves from the
% pose-dependent inflated envelopes (physically simultaneous poses only)
depth = 0;
depthMaxMm = 0;
nB = cs.seed.nBPA;
if nB < 2
    return
end
Drest = fr_diameter(C, cs.diaStr, 0);
exitRow = cs.cross(1);
nG = mech.nG;
for i = 1:nG
    T1 = meta.T1list{i};
    T2 = meta.T2list{i};
    p1 = fr_applyT(T1, meta.loc(1, :));
    p2 = fr_applyT(T2 * T1, meta.loc(exitRow, :));
    u = (p2 - p1) / max(norm(p2 - p1), 1e-9);
    u = u(:)';   % row 1x3 for the perp basis below
    if any(~isfinite(u))
        u = [1; 0; 0]';
    end
    uc = u(:);   % 3x1 column
    perp = null(uc');
    perp = perp(:, 1:min(2, size(perp, 2)));   % 3x1 or 3x2 columns
    if size(perp, 2) < 2 || any(~isfinite(perp(:)))
        % degenerate axis: fixed perpendicular frame, columns 3x1
        alt = [1; 0; 0];
        if abs(dot(alt, uc)) > 0.9
            alt = [0; 1; 0];
        end
        perp = [cross(alt, uc), cross(uc, alt)];   % 3x2, both columns perp
    end
    offs = fr_ringOffsets(nB) * (Drest / 2) / 1000;
    Dinfl = fr_diameter(C, cs.diaStr, mech.strain_p(i));
    r = (Dinfl / 2) / 1000;
    ctr = cell(1, size(offs, 1));
    for b = 1:size(offs, 1)
        o = perp * offs(b, :)';
        ctr{b} = [p1 + o'; p2 + o'];   % 2x3 endpoints
    end
    for a = 1:numel(ctr)
        for b2 = a+1:numel(ctr)
            d = segSegDistMin(ctr{a}(1, :), ctr{a}(2, :), ...
                ctr{b2}(1, :), ctr{b2}(2, :));
            ov = 2 * r - d;
            if ov > 0
                depth = depth + ov / (Drest / 1000);
                depthMaxMm = max(depthMaxMm, ov * 1000);
            end
        end
    end
end
end

function d = segSegDistMin(p1, p2, q1, q2)
u = p2 - p1; v = q2 - q1; w0 = p1 - q1;
a = dot(u, u); b = dot(u, v); c = dot(v, v);
d0 = dot(u, w0); e = dot(v, w0);
D = a * c - b * b;
if D < 1e-24
    sc = 0;
    tc = max(0, min(1, e / max(c, eps)));
else
    sc = max(0, min(1, (b * e - c * d0) / D));
    tc = max(0, min(1, (a * e - b * d0) / D));
    den = a * c - b * b;
    if abs(den) > 1e-24
        sc = max(0, min(1, (b * tc - d0) / den));
        tc = max(0, min(1, (b * sc - e) / max(c, eps)));
    end
end
d = norm(p1 + sc * u - (q1 + tc * v));
end

function offs = fr_ringOffsets(nB)
switch min(nB, 3)
    case 1
        offs = [0 0];
    case 2
        offs = [0 1; 0 -1];
    otherwise
        th = (0:2) * 2 * pi / 3;
        offs = [cos(th)', sin(th)'];
end
end

function p = fr_applyT(Tm, prow)
p = (Tm(1:3, 1:3) * prow(:) + Tm(1:3, 4))';
end

function z = fr_rowZPelvis(cs, kin, q, passIdx, usedRowIdx, locRow)
% express a used-row point (assigned frame) into the pelvis frame at q
frame = cs.rowFrame{passIdx}{fr_absRow(cs, passIdx, usedRowIdx)};
T = fr_chainT(kin, q, frame, 'pelvis');
z = (T(1:3, 1:3) * locRow(:) + T(1:3, 4));
z = z(3);
end

function rAbs = fr_absRow(cs, passIdx, usedRowIdx)
absRows = find(cs.rowUse{passIdx});
rAbs = absRows(usedRowIdx);
end

function margin = fr_fittingTipMargin(C, cs, kin, meta)
% sleeve fitting tip (exit row + fitting length toward the crossing) must
% stay proximal of the ankle pivot plane in the tibia frame (+5 mm).
% The exit row is the last PROXIMAL row, so its assigned frame for the
% foot routes IS the tibia - the tip direction points from the exit
% toward the (transformed) crossing.
ankleY = -0.43;   % tibia_r_offset translation y (osim ankle_r)
exIdx = cs.cross(1) - 1;
margin = inf;
for p = 1:numel(cs.passes)
    qL = cs.passes{p}.qStructs;
    for i = 1:numel(qL)
        T1 = cs.passes{p}.T1(kin, qL{i});
        T2 = cs.passes{p}.T2(kin, qL{i});
        pExit = meta.loc(exIdx, :);            % proximal (tibia) frame
        pCross = fr_applyT(T2 * T1, meta.loc(cs.cross(1), :));
        d = pCross - pExit;
        d = d / max(norm(d), 1e-9);
        tip = pExit + d * C.fittingLength;
        margin = min(margin, tip(2) - ankleY - 0.005);
    end
end
end

function W = fr_rowWorld(cs, kin, q, passIdx, locUsed)
W = zeros(size(locUsed, 1), 3);
mask = cs.rowUse{passIdx};
absRows = find(mask);
T1 = cs.passes{passIdx}.T1(kin, q);
T2 = cs.passes{passIdx}.T2(kin, q);
crossV = cs.crossFor{passIdx};
for u = 1:size(locUsed, 1)
    frame = cs.rowFrame{passIdx}{absRows(u)};
    if u < crossV(1)
        pl = locUsed(u, :);
    elseif u < crossV(2)
        pl = fr_applyT(T1, locUsed(u, :));
    else
        pl = fr_applyT(T2 * T1, locUsed(u, :));
    end
    Tb = fr_bodyT(kin, q, frame);
    W(u, :) = (Tb(1:3, 1:3) * pl(:) + Tb(1:3, 4))';
end
end

function rep = fr_reportShell(rep)
if ~isfield(rep, 'spanRequired')
    rep.spanRequired = [];
end
end

function tt = interpTarget(tt, n)
% resample a target curve to the evaluation grid length (linear)
xi = linspace(0, 1, n);
xo = linspace(0, 1, numel(tt));
tt = interp1(xo, tt, xi, 'linear', 'extrap');
end
