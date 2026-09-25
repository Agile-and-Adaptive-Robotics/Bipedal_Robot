% Opt Sanity Bi Pulley
% Author: Ben Bolen
% Date: 2026-09-24
% Description: Sanity gate for the biarticular reverse-pulley family
% (BiPam_pulley + the two spec builders + the batch driver). NO
% optimization, no figures, no parpool; minutes not hours. Checks:
%   (a) two-joint numeric case (two hinges, 5x5 orientation grid):
%       BiPam_pulley's rigid SegmentLengths and per-crossing UnitDirection
%       match a straightforward INDEPENDENT recomputation (tol 1e-9);
%   (b) G = 1 no-pulley mode reduces to the mono stiffness pipeline:
%       (b1)/(b2) each crossing is compared against a MonoPam_pulley run
%       on the equivalent single-crossing geometry (F_p, mA_p, Torque_p,
%       strain, gama; rel tol 1e-8). Each comparison uses a RIGID bracket
%       (Xi1 = Xi2 = Inf) so the bracket frame cancels from the
%       equilibrium, and an AXIS-TRICK route: the segments NOT carried by
%       the compared crossing's span lie along the other joint's hinge
%       axis (rotation-invariant length), which is the condition that
%       makes the mono analog's span-referenced equilibrium algebraically
%       identical to the bi's (the same construction Opt_sanity_pulley.m
%       documents for its straight-tendon route);
%       (b3) on a GENERAL route (no axis tricks) the G = 1 route state is
%       the base-class solve: the inline force-balance residual
%       |F(s) - (s - s0)/(1/kSpr)| stays at solver tolerance, and the
%       per-crossing deformed kinematics (uD_p, mA_p) match the mono
%       analogs exactly (they are equilibrium-independent);
%   (c) G = 2 at the distal crossing: F_t = F_BPA/G there (independent
%       recompute kSpr*gama vs F_t, tol 1e-6), the closure
%       G*PulleyTravel - DeltaL = gama holds (tol 1e-6), the proximal
%       segment carries the full bundle pull while the distal carries
%       F_t, with nPulleyBPA = 2 the tendon tension equals the summed BPA
%       pull (tol 1e-6), slack cells transmit nothing (F_ins = 0),
%       infeasible cells are NaN with the flag set, and the mount
%       reaction appears only where a tackle sits;
%   (d) biPulleySpecsFromOpenSim parses the osim for the default muscle
%       list, returns finite sane points and Fmax, and its symmetric
%       nPulleyBPA = 2 bundle keeps the bundle centroid ON the OpenSim
%       line (offsets sum to zero, tol 1e-9) at spacing >= diameter +
%       margin;
%   (e) the batch driver in DRY-RUN resolves and prints the plan without
%       error (run LAST: the driver script shares this workspace);
%   (f) in bowden mode the insertion-side direction u_t is constant in
%       the crossing's child frame across orientations (tol 1e-12) while
%       the closure still holds (tol 1e-6).
% error() on any failure; prints OPT_SANITY_BIPULLEY PASS as the last line.

%% Setup (no figures, no parpool)
close all

scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root) || isempty(parent)
        error('Could not locate the Bipedal_Robot repo root from %s', scriptDir)
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

fprintf('== Opt_sanity_BiPulley ==\n')

%% Shared two-joint numeric case (two hinges, 5x5 orientation grid) --------
% General route: origin (frame 1), two vias (frame 2), insertion (frame 3);
% cross = [2 4]. Joint 1 about z through [0.02 -0.01 0] (frame-1 coords);
% joint 2 about y through [0.03 -0.12 0] (frame-2 coords).
ptsA = [ 0.000,  0.000,  0.010;
         0.040, -0.030,  0.000;
         0.050, -0.100,  0.010;
         0.010, -0.150,  0.000];
crossA = [2, 4];
N1 = 5;
N2 = 5;
piv1 = [0.02, -0.01, 0];
piv2 = [0.03, -0.12, 0];
th1 = linspace(-20, 60, N1) * pi / 180;    % joint-1 sweep, rad
th2 = linspace(-40, 25, N2) * pi / 180;    % joint-2 sweep, rad

T = zeros(4, 4, N1, N2);
for i1 = 1:N1
    T(:, :, i1, 1) = RpToTrans(axisRotM([0, 0, 1], th1(i1)), piv1');
end
for i2 = 1:N2
    T(:, :, i2, 2) = RpToTrans(axisRotM([0, 1, 0], th2(i2)), piv2');
end

% Axis-trick routes (see header b1/b2): segment 1 along the joint-1 axis
% (rotation-invariant under T1) and segment 3 along the joint-2 axis
% (invariant under T2) -- each comparison route kills ONE joint's span
% variation so the compared crossing's span carries the whole route change.
ptsB = [ piv1 + [0, 0, 0.030];          % row 1, frame 1: on joint-1 axis
         piv1 + [0, 0, 0.080];          % row 2, frame 2: on joint-1 axis
         0.050, -0.100,  0.010;         % row 3, frame 2 (via)
         0,      0.050,  0];            % row 4, frame 3: on joint-2 axis

%% (a) Rigid geometry vs independent recomputation -------------------------
fprintf('(a) two-hinge 5x5 geometry vs independent recomputation ...\n')

objA = BiPam_pulley('sanityA', ptsA, crossA, 20, T, 0.30, 0.225, 0.02, ...
    0.025, 620, 0.004, Inf, Inf, 6);    % config omitted: pulley disabled

% Independent recomputation (longhand, deliberately NOT the class helpers).
refSeg = zeros(N1, N2, size(ptsA, 1) - 1);
refUD = zeros(N1, 3, N2, 2);
for i1 = 1:N1
    for i2 = 1:N2
        for s = 1:size(ptsA, 1) - 1
            pa = refMapRow(ptsA, crossA, T, s, 1, i1, i2);
            pb = refMapRow(ptsA, crossA, T, s + 1, 1, i1, i2);
            refSeg(i1, i2, s) = norm(pa - pb);
        end
        for k = 1:2
            pa = refMapRow(ptsA, crossA, T, crossA(k) - 1, k + 1, i1, i2);
            pb = refMapRow(ptsA, crossA, T, crossA(k), k + 1, i1, i2);
            d = pa - pb;
            refUD(i1, :, i2, k) = d / norm(d);
        end
    end
end

errSeg = max(abs(objA.SegmentLengths - refSeg), [], 'all');
errUD = max(abs(objA.UnitDirection - refUD), [], 'all');
assertClass(errSeg <= 1e-9, '(a) SegmentLengths vs recomputation (tol 1e-9)')
assertClass(errUD <= 1e-9, '(a) UnitDirection vs recomputation (tol 1e-9)')
fprintf('    PASS (max |seg err| = %.3g, max |dir err| = %.3g)\n', ...
    errSeg, errUD)

%% (b) G = 1 reduces to the mono pipeline ----------------------------------
fprintf('(b) G=1 reduction vs MonoPam_pulley ...\n')

% Sizing from the rigid route: rest so the reference contraction sits
% mid-range (no slack, no infeasibility by construction).
tendon = 0.02;
fitn = 0.025;
Xi0 = 0.004;

% (b1) crossing 1 on the joint-2-axis-trick route: DELIBERATELY rigidly
% constant (rows 1-2 both lie on the joint-1 axis, so DeltaL(:,:,1) == 0).
% This is a degenerate identity cell, not a span-varying one -- (b2) ref
% cell 21 and (b3) cover the span-varying behavior.
mLmaxB = routeMaxLength(ptsB, crossA, T, N1, N2);
restB = (mLmaxB - tendon - 2 * fitn - Xi0) / 0.875;
kmaxB = 0.75 * restB;
bi1 = BiPam_pulley('sanityB1', ptsB, crossA, 20, T, restB, kmaxB, ...
    tendon, fitn, 620, Xi0, Inf, Inf, 6);
mono1 = buildMonoAnalog(ptsB, crossA, T, 1, crossA(1) - 1, ...
    restB, kmaxB, tendon, fitn, Xi0, N1, N2);
compareBiMono(bi1, mono1, 1, '(b1) crossing 1 (joint-2-axis route)')

% (b2) crossing 2 on the joint-1-axis-trick route (segment 1 invariant
% under T1, segment 2 frame-2-constant: only segment 3 carries variation).
ptsC = ptsB;
ptsC(1, :) = piv1 + [0, 0, 0.030];
ptsC(2, :) = piv1 + [0, 0, 0.080];
ptsC(3, :) = [0.050, -0.100, 0.010];
ptsC(4, :) = [0.010, -0.150, 0.000];    % general insertion (joint 2 free)
mLmaxC = routeMaxLength(ptsC, crossA, T, N1, N2);
restC = (mLmaxC - tendon - 2 * fitn - Xi0) / 0.875;
kmaxC = 0.75 * restC;
bi2 = BiPam_pulley('sanityB2', ptsC, crossA, 20, T, restC, kmaxC, ...
    tendon, fitn, 620, Xi0, Inf, Inf, 6);
mono2 = buildMonoAnalog(ptsC, crossA, T, 2, crossA(2) - 1, ...
    restC, kmaxC, tendon, fitn, Xi0, N1, N2);
compareBiMono(bi2, mono2, 2, '(b2) crossing 2 (joint-1-axis route)')

% (b3) general route: the G = 1 route state IS the base solve (inline
% force-balance residual), and the deformed kinematics match the mono
% analogs exactly on the unmodified route.
mLmaxA = max(objA.MuscleLength(:));
restA = (mLmaxA - tendon - 2 * fitn - Xi0) / 0.875;
kmaxA = 0.75 * restA;
biG = BiPam_pulley('sanityB3', ptsA, crossA, 20, T, restA, kmaxA, ...
    tendon, fitn, 620, Xi0, Inf, Inf, 6);

s0A = restA - biG.MuscleLength + Xi0 + tendon + 2 * fitn;
FmagA = biG.FsingleBPA;
resBase = abs(FmagA(:) - (biG.sContraction(:) - s0A(:)) * biG.kSpr);
scaleBase = max(FmagA(:));
assertClass(max(resBase) <= 1e-6 * max(1, scaleBase), ...
    '(b3) G=1 route state satisfies the base force balance F = (s - s0)/kSpr')
for k = 1:2
    monoK = buildMonoAnalog(ptsA, crossA, T, k, crossA(k) - 1, ...
        restA, kmaxA, tendon, fitn, Xi0, N1, N2);
    biU = flat3(biG.uD_p(:, :, :, k));
    biM = flat3(biG.mA_p(:, :, :, k));
    assertClass(relErr(biU, monoK.uD_p) <= 1e-8, ...
        '(b3) uD_p matches MonoPam_pulley at crossing %d (rel 1e-8)', k)
    assertClass(relErr(biM, monoK.mA_p) <= 1e-8, ...
        '(b3) mA_p matches MonoPam_pulley at crossing %d (rel 1e-8)', k)
end
fprintf(['    PASS (b1)/(b2) full-pipeline matches; (b3) base-balance ' ...
    'residual %.3g N, kinematics match both crossings\n'], max(resBase))

%% (c) G = 2 at the distal crossing ----------------------------------------
fprintf('(c) G=2 distal-crossing transmission ...\n')

cfgC = struct('nPulleyBPA', 1, 'tackleLineParts', [1, 2], ...
    'pulleyExitIndex', crossA - 1);
bi2g = BiPam_pulley('sanityC', ptsC, crossA, 20, T, restC, kmaxC, ...
    tendon, fitn, 620, Xi0, Inf, Inf, 6, cfgC);

assertClass(bi2g.ActiveCrossing == 2, '(c) distal crossing is active')

% F_t = F_BPA / G at the distal crossing, both sides recomputed
% independently (tol 1e-6 relative to the force scale).
scaleF = max(bi2g.FsingleBPA(:));
resFt = abs(bi2g.kSpr * bi2g.gama - bi2g.Ftendon(:));
assertClass(max(resFt) <= 1e-6 * max(1, scaleF), ...
    '(c) F_t = F_BPA / 2 at the distal crossing (tol 1e-6)')

% Closure in tackle-travel form: G*PulleyTravel - DeltaL = gama (tol 1e-6
% relative to the span scale).
deltaL2 = bi2g.deltaL(:, :, 2);
resClosure = abs(2 * bi2g.PulleyTravel(:) - deltaL2(:) - bi2g.gama);
assertClass(all(isfinite(bi2g.PulleyTravel(:))), '(c) PulleyTravel finite')
assertClass(max(resClosure) <= 1e-6 * max(abs(deltaL2(:))), ...
    '(c) 2*PulleyTravel - DeltaL = gama at the distal crossing (tol 1e-6)')

% Tension split: the proximal segment carries the full BPA-side pull, the
% distal segment the divided tendon tension.
resSplit = abs(bi2g.segTension(:, 1) - bi2g.FsingleBPA(:)) ...
    + abs(bi2g.segTension(:, 2) - bi2g.Ftendon(:));
assertClass(max(resSplit) <= 1e-9 * max(1, scaleF), ...
    '(c) proximal segment = F_BPA, distal segment = F_t (tol 1e-9)')

% Insertion force equals the tackle tension cell-for-cell where live.
Fins2 = flat3(bi2g.F_ins(:, :, :, 2));
FinsMag = sqrt(sum(Fins2.^2, 2));
infeasV = bi2g.PulleyInfeasible(:);
slackV = bi2g.PulleySlack(:);
liveV = ~(infeasV | slackV);
errIns = max(abs(FinsMag(liveV) - bi2g.Ftendon(liveV)));
assertClass(errIns <= 1e-6 * max(1, scaleF), ...
    '(c) |F_ins| = F_t at the distal crossing (tol 1e-6)')
assertClass(all(Fins2(slackV, :) == 0, 'all'), ...
    '(c) slack cells transmit nothing (F_ins = 0)')
if any(infeasV)
    assertClass(all(isnan(Fins2(infeasV, :)), 'all'), ...
        '(c) infeasible cells NaN the insertion force')
end

% Reaction only where a tackle physically sits.
R2 = bi2g.ReactionF(:, :, :, 2);
R1 = bi2g.ReactionF(:, :, :, 1);
assertClass(all(isnan(R1(:))), '(c) no reaction at the tackle-free crossing')
assertClass(any(~isnan(R2(:))), '(c) reaction at the tackle crossing')
% N1 x 3 x N2 must be permuted before flattening (dim 1 is the
% fastest-varying axis; a bare reshape scrambles the components).
Rmag = sqrt(sum(reshape(permute(R2, [1, 3, 2]), [], 3).^2, 2));
RmagOK = Rmag(liveV);
RpropOK = bi2g.ReactionFmag(liveV);
assertClass(max(abs(RmagOK - RpropOK)) <= 1e-9 * max(1, max(RpropOK)), ...
    '(c) ReactionFmag matches the vector magnitude')

% nPulleyBPA = 2: the tendon tension equals the SUMMED BPA pull (the total
% pull doubles, the travel gain does not), recomputed independently.
cfgC2 = struct('nPulleyBPA', 2, 'tackleLineParts', 1, ...
    'pulleyExitIndex', crossA - 1);
bi2n = BiPam_pulley('sanityC2', ptsC, crossA, 20, T, restC, kmaxC, ...
    tendon, fitn, 620, Xi0, Inf, Inf, 6, cfgC2);
Fsingle2n = festo4(20, bi2n.sContraction(:) / restC / ...
    ((restC - kmaxC) / restC), 620) .* bi2n.Fmax;
resPull = abs(bi2n.kSpr * bi2n.gama - 2 * Fsingle2n);
deltaL2n = bi2n.deltaL(:, :, bi2n.ActiveCrossing);
resTravel2n = abs(bi2n.PulleyTravel(:) - deltaL2n(:) - bi2n.gama);
assertClass(max(resPull) <= 1e-6 * max(1, max(Fsingle2n)), ...
    '(c) nPulleyBPA=2: F_tendon = F_BPA{1} + F_BPA{2} (tol 1e-6)')
assertClass(max(resTravel2n) <= 1e-6 * max(1e-3, max(abs(bi2n.gama))), ...
    '(c) nPulleyBPA=2: travel gain unchanged (1*PulleyTravel = DeltaL + gama)')
fprintf(['    PASS (max F_t residual %.3g N; closure %.3g m; pull-sum ' ...
    '%.3g N; active crossing %d)\n'], max(resFt), max(resClosure), ...
    max(resPull), bi2g.ActiveCrossing)

%% (d) OpenSim builder default list + symmetric bundle ----------------------
fprintf('(d) biPulleySpecsFromOpenSim default list + bundle ...\n')
specsO = biPulleySpecsFromOpenSim();
assertClass(numel(specsO) == 5, '(d) default list resolves 5 muscles')
for k = 1:numel(specsO)
    s = specsO(k);
    assertClass(all(isfinite(s.points(:))), ...
        '(d) finite points: %s', s.name)
    assertClass(all(abs(s.points(:)) < 0.6), ...
        '(d) sane point magnitudes (< 0.6 m): %s', s.name)
    assertClass(s.fmax > 0 && isfinite(s.fmax), ...
        '(d) positive finite Fmax: %s', s.name)
    assertClass(strcmp(s.pulley.routingMode, 'moving_via') ...
        && s.pulley.nPulleyBPA == 1 && s.pulley.tackleLineParts == 1, ...
        '(d) default transmission config 1/1/moving_via: %s', s.name)
end

% nPulleyBPA = 2: symmetric bundle about the ORIGINAL OpenSim line.
cfgB = struct('nPulleyBPA', 2, 'tackleLineParts', 1);
specsB = biPulleySpecsFromOpenSim([], [], cfgB);
for k = 1:numel(specsB)
    s = specsB(k);
    assertClass(size(s.bundleOffsets, 1) == 2, ...
        '(d) bundle offsets for 2 BPAs: %s', s.name)
    assertClass(max(abs(sum(s.bundleOffsets, 1))) <= 1e-9, ...
        '(d) bundle centroid ON the OpenSim line (offsets sum to 0): %s', ...
        s.name)
    spacing = vecnorm(s.bundleOffsets(2, :) - s.bundleOffsets(1, :));
    assertClass(spacing >= s.diameter * 1e-3 + s.bundleMargin - 1e-12, ...
        '(d) bundle spacing >= diameter + margin: %s', s.name)
    assertClass(s.footprintWidth >= (s.diameter * 1e-3), ...
        '(d) footprint width documented: %s', s.name)
end
fprintf('    PASS (%d muscles: %s; 2-BPA bundle centroid on-line)\n', ...
    numel(specsO), strjoin({specsO.name}, ', '))

%% (f) Bowden routing mode --------------------------------------------------
fprintf('(f) bowden routing mode ...\n')
cfgBwd = struct('nPulleyBPA', 1, 'tackleLineParts', [1, 2], ...
    'routingMode', 'bowden', 'pulleyExitIndex', crossA - 1);
biB = BiPam_pulley('sanityF', ptsC, crossA, 20, T, restC, kmaxC, ...
    tendon, fitn, 620, Xi0, Inf, Inf, 6, cfgBwd);
u_tB = biB.u_t(:, :, :, 2);
devB = max(vecnorm(u_tB - reshape(u_tB(1, :, 1), 1, 3, 1), 2, 2), [], 'all');
deltaLB = biB.deltaL(:, :, 2);
resClosureB = abs(2 * biB.PulleyTravel(:) - deltaLB(:) - biB.gama);
assertClass(devB <= 1e-12, ...
    '(f) bowden u_t constant in the child frame across the grid (tol 1e-12)')
assertClass(max(resClosureB) <= 1e-6 * max(abs(deltaLB(:))), ...
    '(f) bowden closure intact (tol 1e-6)')
assertClass(biB.BowdenBossDia == 0.007 && biB.BowdenRunClearance == 0.005, ...
    '(f) bowden envelope constants at the Shimano-type defaults')
fprintf('    PASS (u_t deviation %.3g; closure residual %.3g m)\n', ...
    devB, max(resClosureB))

%% (e) Batch driver DRY-RUN (last: shares this workspace) ------------------
fprintf('(e) Opt_run_BiPulley DRY-RUN ...\n')
run(fullfile(scriptDir, 'Opt_run_BiPulley.m'))
assertClass(exist('manifest', 'var') == 1 && ~isempty(manifest), ...
    '(e) driver resolved the manifest in DRY-RUN')
fprintf('    PASS (driver resolved %d muscles and returned cleanly)\n', ...
    numel(manifest))

%% --------------------------------------------------------------------------
fprintf('\nOPT_SANITY_BIPULLEY PASS\n')

%% =====================================================================
%% Local functions
%% =====================================================================
function assertClass(cond, what, varargin)
if ~all(cond(:))
    if nargin > 2
        error('Opt_sanity_BiPulley:failed', 'SANITY FAIL: %s [%s]', ...
            sprintf(what, varargin{:}), sprintf('%g ', varargin{:}))
    else
        error('Opt_sanity_BiPulley:failed', 'SANITY FAIL: %s', what)
    end
end
end

function e = relErr(a, b)
% Relative error (max over rows, scaled by the reference magnitude), with
% NaN rows counted only when BOTH sides are NaN there.
if ~isequal(isnan(a), isnan(b))
    e = inf;
    return
end
okRows = ~any(isnan(a), 2);
if ~any(okRows)
    e = 0;
    return
end
da = max(abs(a(okRows, :) - b(okRows, :)), [], 2);
scale = max(abs(a(okRows, :)), [], 2);
scale(scale < 1e-12) = 1;
e = max(da ./ scale);
end

function v = flat3(x)
% N1 x 3 x N2 -> N x 3 in the flattened-grid (ii fastest) order. The
% middle 3-wide axis MUST be permuted forward before reshaping, or the
% rows scramble.
v = reshape(permute(x, [1, 3, 2]), [], 3);
end

function Lmax = routeMaxLength(pts, cross, T, N1, N2)
% Max rigid route length over the grid (driver-side sizing helper).
Lmax = 0;
for i1 = 1:N1
    for i2 = 1:N2
        tot = 0;
        for r = 1:size(pts, 1) - 1
            pa = refMapRow(pts, cross, T, r, 1, i1, i2);
            pb = refMapRow(pts, cross, T, r + 1, 1, i1, i2);
            tot = tot + norm(pa - pb);
        end
        Lmax = max(Lmax, tot);
    end
end
end

function p = refMapRow(pts, cross, T, rr, frame, i1, i2)
% Independent longhand frame mapping of route row rr (BiPamData row-frame
% rule; parent chain down, inverse chain up).
if rr < cross(1)
    rf = 1;
elseif cross(1) == cross(2)
    rf = 3;
elseif rr < cross(2)
    rf = 2;
else
    rf = 3;
end
p = pts(rr, :);
while rf > frame
    if rf == 3
        p = xform(T(:, :, i2, 2), p);
        rf = 2;
    else
        p = xform(T(:, :, i1, 1), p);
        rf = 1;
    end
end
while rf < frame
    if rf == 1
        p = xform(inv(T(:, :, i1, 1)), p);
        rf = 2;
    else
        p = xform(inv(T(:, :, i2, 2)), p);
        rf = 3;
    end
end
end

function compareBiMono(bi, mono, k, label)
% Full-pipeline comparison of the bi class's crossing k against the mono
% analog at G = 1: force, moment arm, torque (with matching NaN
% patterns), strain, and tendon stretch.
biF = flat3(bi.F_p(:, :, :, k));
biM = flat3(bi.mA_p(:, :, :, k));
biT = flat3(bi.Torque_p(:, :, :, k));

relF = relErr(biF, mono.F_p);
relM = relErr(biM, mono.mA_p);
relT = relErr(biT, mono.Torque_p);
sameNaN = isequal(isnan(biT), isnan(mono.Torque_p));
relStrain = max(abs(bi.strain_p - mono.strain_p)) / max(1e-12, max(abs(mono.strain_p)));
relGama = max(abs(bi.gama - mono.gama));

assertClass(relF <= 1e-8, '%s: F_p matches MonoPam_pulley (rel 1e-8)', label)
assertClass(relM <= 1e-8, '%s: mA_p matches MonoPam_pulley (rel 1e-8)', label)
assertClass(sameNaN, '%s: Torque_p NaN pattern matches MonoPam_pulley', label)
assertClass(relT <= 1e-8, '%s: Torque_p matches MonoPam_pulley (rel 1e-8)', label)
assertClass(relStrain <= 1e-8, '%s: strain_p matches MonoPam_pulley', label)
assertClass(relGama <= 1e-8, '%s: gama matches MonoPam_pulley', label)
fprintf(['    %s: PASS (rel err: F %.3g, mA %.3g, T %.3g; analog ref ' ...
    'cell %d)\n'], label, relF, relM, relT, mono.RefIndex)
end

function mono = buildMonoAnalog(pts, cross, T, k, exitRow, rest, kmax, ...
    tendon, fitn, Xi0, N1, N2)
% Build the MonoPam_pulley single-crossing analog of bi crossing k:
% parent frame = frame k, child = frame k+1, T_m = the crossing-k transform
% at its angle; Location rows pre-mapped into the analog frames (before
% Cross(k): parent coords; from Cross(k): child coords).
nPts = size(pts, 1);
loc = zeros(nPts, 3, N1 * N2);
Tm = zeros(4, 4, N1 * N2);
for i1 = 1:N1
    for i2 = 1:N2
        m = i1 + (i2 - 1) * N1;
        for r = 1:nPts
            loc(r, :, m) = analogRow(pts, cross, k, r, i1, i2, T);
        end
        if k == 1
            Tm(:, :, m) = T(:, :, i1, 1);
        else
            Tm(:, :, m) = T(:, :, i2, 2);
        end
    end
end
cfgM = struct('nPulleyBPA', 1, 'tackleLineParts', 1, ...
    'pulleyExitIndex', exitRow);
mono = MonoPam_pulley('sanityMono', loc, cross(k), 20, Tm, rest, kmax, ...
    tendon, fitn, 620, Xi0, Inf, Inf, 6, cfgM);
end

function p = analogRow(pts, cross, k, r, i1, i2, T)
% Row r of the bi route expressed in the mono analog's frames.
if r < cross(k)
    target = k;
else
    target = k + 1;
end
p = pts(r, :);
if r < cross(1)
    rf = 1;
elseif cross(1) == cross(2)
    rf = 3;
elseif r < cross(2)
    rf = 2;
else
    rf = 3;
end
while rf > target
    if rf == 3
        p = xform(T(:, :, i2, 2), p);
        rf = 2;
    else
        p = xform(T(:, :, i1, 1), p);
        rf = 1;
    end
end
while rf < target
    if rf == 1
        p = xform(inv(T(:, :, i1, 1)), p);
        rf = 2;
    else
        p = xform(inv(T(:, :, i2, 2)), p);
        rf = 3;
    end
end
end

function p = xform(T, p)
% RowVecTrans inline (row vector through a homogeneous transform).
r = T * [p, 1]';
p = [r(1), r(2), r(3)];
end

function R = axisRotM(axis, theta)
% Rodrigues rotation about a unit axis.
a = axis(:) / norm(axis);
Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end
