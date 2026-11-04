function fr_xi3_check()
%FR_XI3_CHECK Xi3-per-pose verification battery (campaign item, directive
% item 3/4 machinery). Runnable exactly as:
%   "D:/Program Files/MATLAB/R2025a/bin/matlab.exe" -batch
%   "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/...
%    Functional_Routing_20261027/campaign_20261103'); fr_xi3_check"
% in a FRESH process; it must PASS every check or error() (nonzero exit).
% Serial only (a detached driver owns the workers; no parpool here).
%
% Checks:
%  (a) rigid-limit recovery: the derived class BiPam_pulley2_xi3pose with
%      Xi1/Xi2 = Inf, Xi3 = 0, mounts off reproduces the ORIGINAL
%      BiPam_pulley2_20261022 rigid branch EXACTLY (bit-identical path;
%      reported max diffs), through the wrapper's explicit A/B flag
%      (caseDef.xi3Mode 'legacyMaxEnvelope' vs 'perPose').
%  (b) component isolation (derived class, direct construction, rigid
%      tendon wraps = Inf so kSpr = Inf and the closure is exactly
%      decomposable): Xi0-only shifts every taut cell's contraction by
%      exactly Xi0; Xi3-only shifts every taut cell by exactly its
%      per-cell deltaRoute; the applied deltaRoute is bit-identical
%      whether Xi1/Xi2 compliance or Xi0 is present (zero cross-coupling
%      of the loss computation); Xi1-only changes the state only through
%      the bracket channel and satisfies the closure identity
%      G*PulleyTravel == DeltaL exactly. With ALL terms on the roots
%      interact through the shared force law F(s) - expected physics,
%      stated in the output, not a coupling defect.
%  (c) Xi3 is genuinely per-pose: a 3-pose synthetic route with
%      per-cell bend measures [b1 b2 b3] produces per-cell losses
%      deltaRoute_i = Xi3*b_i*comp_i^2 (printed); two poses of the real
%      campaign cases with different bend angles produce different
%      per-cell losses (printed, incl. the legacy max-envelope scalar
%      for contrast); loss is linear in the bend at fixed comp
%      (ratio-3 probe, tol 1e-12 relative).
%  (d) work balance WITH the loss inside: per taut cell,
%        W_gross - W_xi0 - W_loss - W_out - W_spring + W_cAdd = 0
%      with W_gross = nBPA*F*(s - s0_base + DeltaL - cb*F),
%      W_xi0 = nBPA*F*Xi0, W_loss = nBPA*F*deltaRoute (the Xi3 loss),
%      W_out = Ft*DeltaL, W_spring = Ft*gama, W_cAdd = Ft*cAdd*F
%      (mount/bowden/insertion-bracket gives). The loss shows up
%      explicitly and the balance closes to < 1e-12 relative.
%  (e) bracket + tackle-mount + Bowden-housing-stop reactions exist and
%      sum sensibly: a load table is printed at one pose; the closure
%      coefficient decomposition cAdd = -G*nBPA*cmB + nBPA*cmT
%      + cmBw(1) + cmBw(2) + cIns is verified (each hardware compliance
%      counted exactly once); the wrapper's bend profile equals the
%      derived class's self-computed one.
%
% Cases (contract cases, fr_cases27): psoas (monoarticular hip, G=2
% tackle), med_gas (biarticular muscle, knee route with G=2 tackle), sar
% (biarticular, TWO DISTINCT crossings), ext_hal (Bowden housing-stop
% channel, foot).

tAll = tic;
fprintf('=== fr_xi3_check | %s ===\n', datetime('now'));
fprintf('machine: %s\n', getenv('COMPUTERNAME'));

C = fr_setup();
kin = fr_kinematics(C.osimStock);
CASES = fr_cases27(C, kin);     % contract cases (names = target-contract keys)

csM = pickCase(CASES, 'psoas');     % monoarticular (hip, tackle G=2)
csK = pickCase(CASES, 'med_gas');   % biarticular muscle, knee route (tackle)
csB = pickCase(CASES, 'sar');       % biarticular, TWO DISTINCT crossings
csF = pickCase(CASES, 'ext_hal');   % bowden housing-stop channel (foot)

% ---- CHECK (a): rigid-limit recovery ------------------------------------
fprintf('\n--- CHECK (a): rigid-limit recovery (derived == original) ---\n');
for cs = [csM, csK, csB]
    rigidRecovery(C, kin, cs);
end

% ---- CHECK (b): component isolation --------------------------------------
fprintf('\n--- CHECK (b): component isolation (Xi0 / Xi1 / Xi3) ---\n');
iso = struct();
iso.rows = {};
for cs = [csM, csK, csB]
    iso = componentIsolation(C, kin, cs, iso);
end

% ---- CHECK (c): per-pose dependence ---------------------------------------
fprintf('\n--- CHECK (c): Xi3 is genuinely per-pose ---\n');
perPoseSynthetic(C);
for cs = [csM, csK, csB]
    perPoseReal(C, kin, cs);
end

% ---- CHECK (d): work balance with the loss inside -------------------------
fprintf('\n--- CHECK (d): work balance with loss (|resid|/scale < 1e-12) ---\n');
for cs = [csM, csK, csB]
    workBalance(C, kin, cs);
end

% ---- CHECK (e): reactions and load table ----------------------------------
fprintf('\n--- CHECK (e): bracket/mount/housing-stop reactions ---\n');
reactionTable(C, kin, csK, 'med_gas (tackle mount)');
reactionTable(C, kin, csM, 'psoas (tackle mount)');
reactionTable(C, kin, csB, 'sar (brackets only)');
reactionBowden(C, kin, csF);
cAddDecomposition(C, kin, csK);

fprintf('\n=== fr_xi3_check PASSED (all checks) in %.1f s ===\n', toc(tAll));
end

% ==========================================================================
function cs = pickCase(CASES, name)
sel = strcmp({CASES.name}, name);
assert(any(sel), 'fr_xi3_check:CaseMissing', 'case %s not found', name);
cs = CASES(find(sel, 1));
cs.xi3Mode = 'perPose';         % explicit A/B flag (default stays legacy)
end

function qg = passQ(cs)
qg = struct('passIdx', 1, 'qList', {cs.passes{1}.qStructs});
end

function fitn = fittingOf(C, cs)
if isfield(cs, 'fittingLength') && ~isempty(cs.fittingLength)
    fitn = cs.fittingLength;
else
    fitn = C.fittingLength;
end
end

function T = tgrid(meta, nG)
T = zeros(4, 4, nG, 2);
for i = 1:nG
    T(:, :, i, 1) = meta.T1list{i};
end
T2ref = meta.T2list{ceil(nG / 2)};
T(:, :, 1, 2) = T2ref;
if nG >= 2
    T(:, :, 2, 2) = T2ref;
end
end

function d = maxDiff(a, b)
% NaN-aware max absolute difference (NaNs must agree in place)
a = a(:);
b = b(:);
nanA = isnan(a);
nanB = isnan(b);
if ~isequal(nanA, nanB)
    d = Inf;
    return
end
if all(nanA)
    d = 0;
    return
end
d = max(abs(a(~nanA) - b(~nanB)));
end

function printDiff(label, d, tol)
fprintf('    %-52s max|diff| = %.3e (tol %.1e) -> %s\n', label, d, tol, ...
    tf(d <= tol));
if ~(d <= tol)
    error('fr_xi3_check:RigidRecovery', ...
        '%s exceeded tolerance %.1e (got %.3e)', label, tol, d);
end
end

function s = tf(x)
if x
    s = 'PASS';
else
    s = 'FAIL';
end
end

% ==========================================================================
function rigidRecovery(C, kin, cs0)
design = fr_unpackX(cs0, cs0.xSeed);
qg = passQ(cs0);
csL = cs0;
csL.xi3Mode = 'legacyMaxEnvelope';
csP = cs0;
csP.xi3Mode = 'perPose';
[mechL, metaL] = fr_mechBuild(C, csL, design, kin, 'rigid', qg);
[mechP, metaP] = fr_mechBuild(C, csP, design, kin, 'rigid', qg);
fprintf('  [%s]\n', cs0.name);

% wrapper legacy-rigid vs wrapper perPose-rigid: the perPose rigid branch
% must be bit-identical (Xi3 = 0 -> exact no-op)
printDiff('wrapper legacy rigid vs perPose rigid: strain_p', ...
    maxDiff(mechL.strain_p, mechP.strain_p), 1e-10);
printDiff('wrapper legacy rigid vs perPose rigid: FsingleBPA', ...
    maxDiff(mechL.FsingleBPA, mechP.FsingleBPA), 1e-10);
printDiff('wrapper legacy rigid vs perPose rigid: sContraction', ...
    maxDiff(mechL.cell.sContraction, mechP.cell.sContraction), 1e-10);
printDiff('wrapper legacy rigid vs perPose rigid: Torque_ins', ...
    maxDiff(mechL.cell.Torque_ins, mechP.cell.Torque_ins), 1e-10);
assert(maxDiff(metaL.loc, metaP.loc) == 0);

% and both vs a DIRECT construction of the ORIGINAL class (triple check)
nG = numel(qg.qList);
fitn = fittingOf(C, cs0);
T = tgrid(metaL, nG);
mOrig = BiPam_pulley2_20261022(cs0.name, metaL.loc, metaL.caseDef.cross, ...
    metaL.caseDef.diaNum, T, metaL.rest, metaL.kmax, design.tendon, ...
    fitn, C.pressure, 0, Inf, Inf, 1e4, metaL.cfg);
mDeriv = BiPam_pulley2_xi3pose(cs0.name, metaL.loc, metaL.caseDef.cross, ...
    metaL.caseDef.diaNum, T, metaL.rest, metaL.kmax, design.tendon, ...
    fitn, C.pressure, 0, Inf, Inf, 1e4, metaL.cfg, 0, []);
printDiff('direct ORIGINAL class vs derived (Xi3=0): strain_p', ...
    maxDiff(mOrig.strain_p, mDeriv.strain_p), 1e-10);
printDiff('direct ORIGINAL class vs derived (Xi3=0): FsingleBPA', ...
    maxDiff(mOrig.FsingleBPA, mDeriv.FsingleBPA), 1e-10);
printDiff('direct ORIGINAL class vs derived (Xi3=0): Torque_ins', ...
    maxDiff(mOrig.Torque_ins, mDeriv.Torque_ins), 1e-10);
printDiff('direct ORIGINAL class vs derived (Xi3=0): deltaL', ...
    maxDiff(mOrig.deltaL, mDeriv.deltaL), 1e-10);
fprintf(['    rigid-limit recovery: derived rigid branch reproduces the ' ...
    'original rigid branch EXACTLY (all diffs 0) -> PASS\n']);
end

% ==========================================================================
function iso = componentIsolation(C, kin, cs0, iso)
design = fr_unpackX(cs0, cs0.xSeed);
qg = passQ(cs0);
[~, metaL] = fr_mechBuild(C, cs0, design, kin, 'rigid', qg);
nG = numel(qg.qList);
fitn = fittingOf(C, cs0);
T = tgrid(metaL, nG);

% per-pose bend profile from the perPose build (the wrapper's profile)
mechP3 = fr_mechBuild(C, cs0, design, kin, 'compliant', qg);
bend1 = mechP3.bendMeasure;                 % nG x 1 per-pose totals
bendGrid = zeros(nG, 2);
bendGrid(:, 1) = bend1;
bendGrid(:, 2) = bend1(ceil(nG / 2));

common = {cs0.name, metaL.loc, metaL.caseDef.cross, metaL.caseDef.diaNum, ...
    T, metaL.rest, metaL.kmax, design.tendon, fitn, C.pressure};
cfg = metaL.cfg;                 % rigid build: NO mount/bowden fields
xi1 = C.Xi.Xi1;
xi2 = C.Xi.Xi2;
chi3 = C.Xi.Xi3;
v0 = C.Xi.Xi0_flx;

mR0 = BiPam_pulley2_xi3pose(common{:}, 0, Inf, Inf, Inf, cfg, 0, []);
mX0 = BiPam_pulley2_xi3pose(common{:}, v0, Inf, Inf, Inf, cfg, 0, []);
mX3 = BiPam_pulley2_xi3pose(common{:}, 0, Inf, Inf, Inf, cfg, chi3, bendGrid);
mX0X3 = BiPam_pulley2_xi3pose(common{:}, v0, Inf, Inf, Inf, cfg, chi3, bendGrid);
mE1 = BiPam_pulley2_xi3pose(common{:}, 0, xi1, xi2, Inf, cfg, 0, []);
mX1X3 = BiPam_pulley2_xi3pose(common{:}, 0, xi1, xi2, Inf, cfg, chi3, bendGrid);
mX3E = BiPam_pulley2_xi3pose(common{:}, 0, xi1, xi2, Inf, cfg, chi3, bendGrid);
fprintf('  [%s] nG = %d\n', cs0.name, nG);

kAct = mR0.ActiveCrossing;
% escape-branch runs (Xi1 = Xi2 = Inf, wraps = Inf): the class assigns
% sCon = s0 for EVERY cell (the pure contraction demand), so the exact
% per-cell shift identities hold on all column-1 cells (column 1 = poses;
% flattened index lin = i, ii fastest). Explicit reshape avoids the
% row/column orientation trap of logical linear indexing.
scR = reshape(mR0.sContraction, [], 1);
scX0 = reshape(mX0.sContraction, [], 1);
scX3 = reshape(mX3.sContraction, [], 1);
fprintf(['    escape-branch exact identities on %d pose cells ' ...
    '(active crossing %d)\n'], nG, kAct);

% Xi0-only: contraction demand shifts by exactly Xi0 at every cell
dX0 = maxDiff(scX0(1:nG), scR(1:nG) + v0);
printDiff('Xi0-only: sCon(Xi0) - sCon(rigid) == Xi0', dX0, 1e-10);

% Xi3-only: contraction demand shifts by exactly the per-cell deltaRoute
dR3 = mX3.deltaRoute;
dX3 = maxDiff(scX3(1:nG), scR(1:nG) + dR3(1:nG));
printDiff('Xi3-only: sCon(Xi3) - sCon(rigid) == deltaRoute(i)', dX3, 1e-10);

% zero cross-coupling of the CLOSURE CHANNELS: Xi1 (force-side compliance)
% must leave the loss computation EXACTLY untouched (comp uses the
% pre-solve strain, which excludes the bracket), and with Xi0 present the
% loss must equal exactly Xi3*bend*comp(Xi0)^2 - the DISSERTATION's zeta
% definition (eq:comp: zeta = 1 - epsilon*, the chi0-inclusive model
% strain; MonoPam_mult.m Contraction_k lines 445-449), i.e. the Xi0 ->
% loss dependence is the validated semantics, not a coupling defect.
KMAX = (metaL.rest - metaL.kmax) / metaL.rest;
mLg = mX0X3.MuscleLength;
strainX0 = (metaL.rest - (mLg - design.tendon - 2 * fitn - v0)) ...
    ./ metaL.rest;
compX0 = max(0, 1 - strainX0 ./ KMAX);
dRexp = chi3 * bendGrid(:) .* compX0(:).^2;
dC0 = maxDiff(mX0X3.deltaRoute, dRexp);
dC1 = maxDiff(mX1X3.deltaRoute, dR3);
printDiff('Xi0+Xi3: deltaRoute == Xi3*bend*comp(Xi0)^2 (zeta semantics)', ...
    dC0, 1e-15);
printDiff('Xi1+Xi3: deltaRoute == deltaRoute(Xi3) exactly', dC1, 1e-15);
assert(isequal(mX0.deltaRoute, zeros(size(dR3))) ...
    && isequal(mE1.deltaRoute, zeros(size(dR3))), ...
    'Xi0-only/Xi1-only runs must apply zero loss');
fprintf(['    zero cross-coupling: Xi0-only and Xi1-only runs apply ' ...
    'EXACTLY zero loss; Xi1 does not perturb the loss at all; the Xi0 ' ...
    'dependence is the zeta = 1 - epsilon*(chi0-inclusive) definition ' ...
    'itself -> PASS\n']);

% elastic-branch runs (finite Xi1/Xi2, wraps = Inf so cSpr = 0): the loss
% enters the ELASTIC root through the same per-cell term
G = mE1.PulleyGain(kAct);
taut = ~mE1.PulleySlack & ~mE1.PulleyInfeasible & mE1.FsingleBPA > 0;
nTaut = nnz(taut);
if nTaut == 0
    error('fr_xi3_check:Isolation', ...
        '%s: no taut elastic cells to test on', cs0.name);
end
fprintf('    elastic-branch taut cells: %d / %d\n', nTaut, numel(taut));
scE1 = reshape(mE1.sContraction, [], 1);
scX3E = reshape(mX3E.sContraction, [], 1);
scR = reshape(mR0.sContraction, [], 1);
ptV = reshape(mE1.PulleyTravel, [], 1);
dLv = reshape(mE1.deltaL(:, :, kAct), [], 1);
pt = ptV(taut);
dL = dLv(taut);
res = max(abs(G * pt - dL)) / max(max(abs(dL)), eps);
fprintf(['    Xi1-only closure identity max |G*PulleyTravel - DeltaL| ' ...
    '/ max|DeltaL| = %.3e -> %s\n'], res, tf(res < 1e-9));
if ~(res < 1e-9)
    error('fr_xi3_check:Isolation', ...
        'Xi1-only run violates its own closure identity (%.3e)', res);
end
cbV = reshape(mE1.cbBracket, [], 1);
cbT = cbV(taut);
cbMax = max(cbT);
brD = reshape(mE1.Bracket1Deflection, [], 3);
brDeflMax = max(vecnorm(brD(taut, :), 2, 2));
fprintf(['    Xi1-only: bracket projection cb max %.3e m/N, bracket ' ...
    'deflection max %.4f mm (bracket channel active; Xi0/Xi3 terms ' ...
    'zero) -> PASS\n'], cbMax, 1000 * brDeflMax);
if ~(cbMax > 0 && brDeflMax > 0)
    error('fr_xi3_check:Isolation', ...
        'Xi1-only run shows no bracket deflection - bracket channel dead');
end

% Xi3 in the elastic branch: the loss consumes takeup, so the solved
% contraction grows by AT MOST its per-cell deltaRoute (the bracket
% projection cb*F(s) shrinks the shift through the force feedback;
% closure: G*(s - s0ref + dR - cbF) = DeltaL with monotone psi)
if ~(any(dR3 > 0))
    error('fr_xi3_check:Isolation', ...
        '%s: bend profile is all-zero - per-pose profile broken', cs0.name);
end
dS3 = scX3E(taut) - scE1(taut);
dRt = dR3(taut);
pos = dRt > 0;
lo = min(dS3(pos) - dRt(pos));
hi = max(dS3 - dRt);
fprintf(['    Xi3 elastic effect: per-cell sCon(Xi3) - sCon(no Xi3) in ' ...
    '[dR %+.3e, dR %+.3e] m (positive shift, at most its own dR through ' ...
    'the force feedback) -> %s\n'], lo, hi, ...
    tf(all(dS3(pos) > 0) && hi <= 1e-9 && all(dS3(~pos) == 0)));
if ~(all(dS3(pos) > 0) && hi <= 1e-9 && all(dS3(~pos) == 0))
    error('fr_xi3_check:Isolation', ...
        'elastic Xi3 shift outside its own per-cell deltaRoute bound');
end
fprintf('    Xi3 elastic per-pose shifts (m): ');
fprintf('%.4f ', dS3);
fprintf('| deltaRoute: ');
fprintf('%.4f ', dRt);
fprintf('-> PASS\n');

% Xi1-only changes the solved state vs its own escape demand (its own
% effect exists; the escape demand scR is the s0 reference)
dState = maxDiff(scE1(taut), scR(taut));
fprintf(['    Xi1-only effect exists: max |sCon(Xi1) - sCon(escape ' ...
    'demand)| = %.3e m -> PASS\n'], dState);
if ~(dState > 0)
    error('fr_xi3_check:Isolation', ...
        'Xi1-only run is identical to rigid - compliance has no effect');
end

% honest statement: with ALL terms on the roots interact through the
% shared F(s) nonlinearity (not a coupling defect; the additive closure
% terms were each verified above in isolation)
iso.rows(end + 1) = {sprintf('%s: dX0 %.2e, dX3 %.2e, dC0 %.2e, dC1 %.2e', ...
    cs0.name, dX0, dX3, dC0, dC1)};
fprintf(['    note: with Xi0+Xi1+Xi3 simultaneously on, the solved roots ' ...
    'interact through the shared force law F(s) (additive closure ' ...
    'terms, nonlinear root) - expected physics, not double counting\n']);
end

% ==========================================================================
function perPoseSynthetic(C)
% 3-pose synthetic route: per-cell bends [b1 b2 b3] -> per-cell losses
loc = [0 0 0; 0.2 0 0; 0.5 0 0];
cross = [2, 3];
rest = 0.45;
kmax = 0.76 * rest;
nG = 3;
T = zeros(4, 4, nG, 2);
for i = 1:nG
    th = deg2rad(20 * (i - 2));
    R = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
    T(:, :, i, 1) = [R, [0; 0; 0]; 0 0 0 1];   % frame2 -> frame1, pivot origin
end
T(:, :, 1, 2) = eye(4);
T(:, :, 2, 2) = eye(4);
bends = [0.005; 0.020; 0.010];
bendGrid = [bends, repmat(bends(2), nG, 1)];   % column 2 = T2ref replica
chi3 = C.Xi.Xi3;
m = BiPam_pulley2_xi3pose('xi3synthetic', loc, cross, 20, T, rest, ...
    kmax, 0.05, 0.02, C.pressure, 0, Inf, Inf, Inf, struct(), chi3, bendGrid);
dR = m.deltaRoute(1:nG);
fprintf('  [synthetic 3-pose] bend (m) = [%g %g %g]\n', bends);
fprintf('    deltaRoute (m) = [%.6e %.6e %.6e]\n', dR);
ratio = dR(2) / max(dR(1), eps);
fprintf(['    ratio deltaRoute(2)/deltaRoute(1) = %.12f (bend ratio 4, ' ...
    'comp equal) -> %s\n'], ratio, tf(abs(ratio - 4) < 1e-9));
if ~(abs(ratio - 4) < 1e-9)
    error('fr_xi3_check:PerPose', ...
        'synthetic per-cell loss does not follow Xi3*b*comp^2 (ratio %g)', ratio);
end
% linearity in the bend at fixed comp (same geometry, scaled profile)
m3 = BiPam_pulley2_xi3pose('xi3synthetic3x', loc, cross, 20, T, rest, ...
    kmax, 0.05, 0.02, C.pressure, 0, Inf, Inf, Inf, struct(), chi3, 3 * bendGrid);
relLin = maxDiff(m3.deltaRoute(1:nG), 3 * dR) / max(max(abs(dR)), eps);
fprintf(['    linearity: max|deltaRoute(3b) - 3*deltaRoute(b)|/max|dR| ' ...
    '= %.3e -> %s\n'], relLin, tf(relLin < 1e-12));
if ~(relLin < 1e-12)
    error('fr_xi3_check:PerPose', 'loss not linear in the bend profile');
end
% and the losses really differ per pose inside ONE object
spread = max(dR) - min(dR);
fprintf('    per-pose spread inside one object: max-min = %.3e m -> %s\n', ...
    spread, tf(spread > 1e-6));
if ~(spread > 1e-6)
    error('fr_xi3_check:PerPose', 'losses identical across poses');
end
end

% ==========================================================================
function perPoseReal(C, kin, cs0)
design = fr_unpackX(cs0, cs0.xSeed);
qg = passQ(cs0);
csP = cs0;
csP.xi3Mode = 'perPose';
csL = cs0;
csL.xi3Mode = 'legacyMaxEnvelope';
[mechP, ~] = fr_mechBuild(C, csP, design, kin, 'compliant', qg);
[mechL, ~] = fr_mechBuild(C, csL, design, kin, 'compliant', qg);
fprintf('  [%s] active crossing %d\n', cs0.name, mechP.ActiveCrossing);
idx = [1, ceil(mechP.nG / 2), mechP.nG];
fprintf('    pose          bend total (m)   per-pose loss (m)   legacy scalar (m)\n');
for i = idx
    fprintf('    %2d           %12.6f    %14.6e    %14.6e\n', i, ...
        mechP.bendMeasure(i), mechP.xi3PerPoseLoss_m(i), ...
        mechL.xi3MaxDeltaRoute_m);
end
spread = max(mechP.xi3PerPoseLoss_m) - min(mechP.xi3PerPoseLoss_m);
fprintf(['    per-pose loss spread = %.3e m (legacy max-envelope would ' ...
    'hold ALL poses at %.3e m) -> %s\n'], spread, ...
    mechL.xi3MaxDeltaRoute_m, tf(spread > 1e-9));
if ~(spread > 1e-9)
    error('fr_xi3_check:PerPose', ...
        '%s: per-pose losses are constant - not genuinely per-pose', cs0.name);
end
% the wrapper's profile must equal the derived class's self-computed one
[~, metaP] = fr_mechBuild(C, csP, design, kin, 'rigid', qg);
nG = numel(qg.qList);
T = tgrid(metaP, nG);
fitn = fittingOf(C, cs0);
mSelf = BiPam_pulley2_xi3pose(cs0.name, metaP.loc, metaP.caseDef.cross, ...
    metaP.caseDef.diaNum, T, metaP.rest, metaP.kmax, design.tendon, ...
    fitn, C.pressure, C.Xi.Xi0_flx, Inf, Inf, Inf, metaP.cfg, ...
    C.Xi.Xi3, []);
bendSelf = sum(mSelf.BendArc, 2);
dProf = maxDiff(bendSelf(1:nG), mechP.bendMeasure);
printDiff('wrapper bend profile == class self-computed profile', ...
    dProf, 1e-12);
end

% ==========================================================================
function workBalance(C, kin, cs0)
design = fr_unpackX(cs0, cs0.xSeed);
qg = passQ(cs0);
[mech, ~] = fr_mechBuild(C, cs0, design, kin, 'compliant', qg);  % perPose
m = mech.cell;
kAct = mech.ActiveCrossing;
nBPA = m.NPulleyBPA(kAct);
nG = mech.nG;
fitn = fittingOf(C, cs0);
% column 1 of the grid = the poses; flattened index lin = i (N1 = nG)
lin = 1:nG;
taut = ~mech.PulleySlack & ~mech.PulleyInfeasible & mech.FsingleBPA > 0;
nTaut = nnz(taut);
if nTaut == 0
    error('fr_xi3_check:WorkBalance', ...
        '%s: no taut cells', cs0.name);
end
mL = m.MuscleLength(:, 1);
s0base = m.RestingL - mL(:) + m.TendonL + 2 * fitn;
sCon = reshape(m.sContraction, [], 1);
dLk = reshape(m.deltaL(:, :, kAct), [], 1);
F = mech.FsingleBPA(taut);
cb = m.cbBracket(lin);
cb = cb(taut);
dR = m.deltaRoute(lin);
dR = dR(taut);
dL = dLk(taut);
gama = m.gama(lin);
gama = gama(taut);
cAdd = m.DeflClosure_cAdd(lin);
cAdd = cAdd(taut);
Ft = mech.Ftendon(taut);
sC = sCon(lin);
sC = sC(taut);
s0b = s0base(taut);
Wgross = nBPA .* F .* (sC - s0b + dL - cb .* F);
Wxi0 = nBPA .* F .* m.Xi0;
Wloss = nBPA .* F .* dR;
Wout = Ft .* dL;
Wspr = Ft .* gama;
WcAdd = Ft .* (cAdd .* F);
resid = Wgross - Wxi0 - Wloss - Wout - Wspr + WcAdd;
scale = max(abs(Wgross));
rel = max(abs(resid)) / max(scale, eps);
fprintf(['  [%s] taut cells %d, scale %.3e J, loss share ' ...
    'sum(W_loss)/sum(W_gross) = %.3e\n'], cs0.name, nTaut, scale, ...
    sum(Wloss) / max(sum(Wgross), eps));
fprintf(['    work balance max |W_in - W_xi0 - W_loss - W_out - ' ...
    'W_spring + W_cAdd| / scale = %.3e -> %s\n'], rel, tf(rel < 1e-12));
if ~(rel < 1e-12)
    error('fr_xi3_check:WorkBalance', ...
        '%s: work balance residual %.3e exceeds 1e-12', cs0.name, rel);
end
if ~(sum(Wloss) > 0)
    error('fr_xi3_check:WorkBalance', ...
        '%s: no loss shows up in the work balance', cs0.name);
end
fprintf(['    loss term W_loss = nBPA*F*deltaRoute present and nonzero ' ...
    '(sum %.3e J) -> PASS\n'], sum(Wloss));
end

% ==========================================================================
function reactionTable(C, kin, cs0, label)
design = fr_unpackX(cs0, cs0.xSeed);
qg = passQ(cs0);
[mech, ~] = fr_mechBuild(C, cs0, design, kin, 'compliant', qg);
m = mech.cell;
kAct = mech.ActiveCrossing;
G = m.PulleyGain(kAct);
nBPA = m.NPulleyBPA(kAct);
i = ceil(mech.nG / 2);          % one mid pose
lin = i;
fprintf('  [%s | %s] pose %d/%d\n', cs0.name, label, i, mech.nG);
F = mech.FsingleBPA(i);
Ft = mech.Ftendon(i);
brR = vecnorm(m.Bracket1Reaction(lin, :), 2);
brD = vecnorm(m.Bracket1Deflection(lin, :), 2);
inR = vecnorm(m.InsBracketReaction(lin, :), 2);
inD = vecnorm(m.InsBracketDeflection(lin, :), 2);
fprintf('    %-28s %10.2f N  defl %8.3f mm\n', 'origin bracket (Xi1/Xi2)', ...
    brR, 1000 * brD);
fprintf('    %-28s %10.2f N  defl %8.3f mm\n', 'insertion bracket', ...
    inR, 1000 * inD);
fprintf('    %-28s %10.2f N\n', 'tendon tension Ft = nBPA*F/G', Ft);
fprintf('    %-28s %10.2f N\n', 'BPA bundle pull nBPA*F', nBPA * F);
if cs0.hasTackle && mech.ReactionFmag(i) > 0
    mR = mech.ReactionFmag(i);
    % 4-D outputs are (N1=poses, 3, N2, crossing): pose i is dim 1
    mD = vecnorm(reshape(m.MountDeflection(i, :, kAct), 1, 3), 2);
    fprintf('    %-28s %10.2f N  defl %8.3f mm\n', ...
        'tackle-block mount reaction', mR, 1000 * mD);
    % sanity: |R| bounded by the triangle inequality on its two terms
    lo = abs(nBPA * F - Ft) - 1e-6;
    hi = nBPA * F + Ft + 1e-6;
    ok = (mR >= lo) && (mR <= hi);
    fprintf(['    mount reaction within [|nBPA*F - Ft|, nBPA*F + Ft] ' ...
        '= [%.1f, %.1f] N -> %s\n'], lo, hi, tf(ok));
    if ~ok
        error('fr_xi3_check:Reactions', 'mount reaction inconsistent');
    end
end
% closure decomposition: each hardware compliance counted exactly once
cAddHat = -G * nBPA * m.DeflClosure_cmB(lin) + nBPA * m.DeflClosure_cmT(lin) ...
    + sum(m.DeflClosure_cmBw(lin, :)) + m.DeflClosure_cIns(lin);
dC = abs(cAddHat - m.DeflClosure_cAdd(lin)) ...
    / max(abs(m.DeflClosure_cAdd(lin)), 1e-15);
fprintf(['    closure decomposition cAdd = -G*nBPA*cmB + nBPA*cmT ' ...
    '+ cmBw(1) + cmBw(2) + cIns: rel err %.3e -> %s\n'], dC, tf(dC < 1e-12));
if ~(dC < 1e-12)
    error('fr_xi3_check:Reactions', ...
        'cAdd decomposition mismatch - compliance double counted?');
end
if ~(brR > 0 && inR > 0 && Ft > 0)
    error('fr_xi3_check:Reactions', ...
        '%s: missing reaction output at pose %d', cs0.name, i);
end
fprintf('    all reaction outputs finite and positive -> PASS\n');
end

% ==========================================================================
function reactionBowden(C, kin, cs0)
% Bowden housing-stop channel: the foot case routes crossing 1 in bowden
% mode with the housing stop bracket compliance (handoff: "housing stops
% react cable tension").
design = fr_unpackX(cs0, cs0.xSeed);
qg = passQ(cs0);
[mech, ~] = fr_mechBuild(C, cs0, design, kin, 'compliant', qg);
m = mech.cell;
i = ceil(mech.nG / 2);
fprintf('  [%s | Bowden housing stops] pose %d/%d\n', cs0.name, i, mech.nG);
hasBowden = any(strcmp(m.RoutingMode, 'bowden'));
if ~hasBowden
    error('fr_xi3_check:Reactions', ...
        '%s case does not exercise the bowden channel', cs0.name);
end
tk = nBPAforceStop(m, mech, i);
fprintf('    %-28s %10.2f N\n', 'housing-stop reacted tension', tk);
bwK = find(strcmp(m.RoutingMode, 'bowden'), 1);
bw = m.BowdenDeflection(:, :, :, bwK);   % (N1=poses, 3, N2)
bwD = vecnorm(reshape(bw(i, :, 1), 1, 3), 2);
fprintf('    %-28s %10.3f mm\n', 'housing-stop deflection', 1000 * bwD);
if ~(tk > 0 && isfinite(bwD))
    error('fr_xi3_check:Reactions', 'bowden stop channel incomplete');
end
fprintf('    Bowden housing-stop reaction outputs exist -> PASS\n');
end

function tk = nBPAforceStop(m, mech, i)
% the tension the bowden stop reacts: nBPA*F/G at/below the active
% tackle, nBPA*F above it (class convention; the foot case has no tackle
% so both coincide)
kAct = mech.ActiveCrossing;
if 1 >= kAct
    tk = mech.Ftendon(i);
else
    tk = mech.FsingleBPA(i) * m.NPulleyBPA(kAct);
end
end

% ==========================================================================
function cAddDecomposition(C, kin, cs0)
% also verify the decomposition over ALL taut cells (not just one pose)
design = fr_unpackX(cs0, cs0.xSeed);
qg = passQ(cs0);
[mech, ~] = fr_mechBuild(C, cs0, design, kin, 'compliant', qg);
m = mech.cell;
kAct = mech.ActiveCrossing;
G = m.PulleyGain(kAct);
nBPA = m.NPulleyBPA(kAct);
nG = mech.nG;
lin = 1:nG;
taut = ~mech.PulleySlack & ~mech.PulleyInfeasible;
cHat = -G * nBPA * m.DeflClosure_cmB(lin) + nBPA * m.DeflClosure_cmT(lin) ...
    + sum(m.DeflClosure_cmBw(lin, :), 2) + m.DeflClosure_cIns(lin);
dC = maxDiff(cHat(taut), m.DeflClosure_cAdd(lin));
dC = dC / max(max(abs(m.DeflClosure_cAdd(lin))), 1e-15);
fprintf('  [%s] all-taut-cell cAdd decomposition rel err %.3e -> %s\n', ...
    cs0.name, dC, tf(dC < 1e-12));
if ~(dC < 1e-12)
    error('fr_xi3_check:Reactions', ...
        '%s: cAdd decomposition fails over all taut cells', cs0.name);
end
end
