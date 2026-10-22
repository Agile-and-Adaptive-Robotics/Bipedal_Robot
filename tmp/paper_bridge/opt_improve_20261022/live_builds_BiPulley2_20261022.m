% live_builds_BiPulley2_20261022 -- one short live build for each new
% mechanics path of the dated variant class BiPam_pulley2_20261022
% (directive PULLEY_SWEEP_DIRECTIVE_20261021.md sections 2 and 3):
%
%   Build A (deflection path): moving-via tackle-block mount deflection
%     from the BPA muscle origin, tackle-increased internal force carried
%     through the equilibrium (mountMaterial aluminum6061T6 + explicit
%     section), plus the aluminum/CroMoly stiffness-ratio check
%     (E 68.9 vs 205 GPa).
%   Build B (bowden path): routingMode 'bowden' with the housing anchor
%     bracket deflection at CroMoly stiffness, pinned-u_t check, and the
%     packaging-envelope fields (boss 7 mm / clearance 5 mm, tunable).
%
% Route/transform source = the recorded doctrine fixture med_gas_r
% (2 BPAs, tackleLineParts [2 1], moving_via), read-only; ONE class
% construction per build -- no optimizer, no parpool, seconds.
% Run: matlab -batch "run('tmp/paper_bridge/opt_improve_20261022/live_builds_BiPulley2_20261022.m')"

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

fixture = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_opensim_med_gas_r_2BPA_G2_20261021_1114.mat');
S = load(fixture);
rest = S.xBest(1);
tendon = S.xBest(2);
kmax = 0.75 * rest;
secStandIn = [0.025, 0.005, 0.040];    % [b h L] m plate, stand-in dims (explicit here)

build = @(cfg) BiPam_pulley2_20261022(S.spec.name, S.points, S.spec.cross, ...
    20, S.T, rest, kmax, tendon, 0.025, S.basePressure, S.baseXi0, ...
    S.baseXi1, S.baseXi2, S.baseWraps, cfg);

% Zero-compliance twin (the default = base-class arithmetic exactly).
o0 = build(S.spec.pulley);

%% ---------------- Build A: moving-via mount deflection ---------------
fprintf('=== Build A: moving-via mount deflection (aluminum 6061-T6) ===\n');
cfgA = S.spec.pulley;
cfgA.mountMaterial = 'aluminum6061T6';
cfgA.mountSection = secStandIn;
oA = build(cfgA);
aK = oA.ActiveCrossing;

assert(any(~isnan(oA.MountDeflection(:))), 'A: MountDeflection all-NaN')
assert(all(isfinite(oA.MountDeflection(~isnan(oA.MountDeflection)))), ...
    'A: MountDeflection not finite')
assert(any(oA.DeflClosure_cAdd ~= 0), 'A: closure coefficient cAdd identically 0')
dA = maxDefl(oA.MountDeflection(:, :, :, aK));
fprintf('  active crossing = %d, tackle nBPA = %d, G = %g\n', aK, ...
    oA.NPulleyBPA(aK), oA.PulleyGain(aK));
fprintf('  max |mount deflection| = %.4f mm; max ReactionFmag = %.1f N\n', ...
    1000 * dA, max(oA.ReactionFmag(:)));
m0 = minTau(o0, aK);
mA = minTau(oA, aK);
dT = maxDeltaTau(o0, oA, aK);
fprintf(['  min |Torque_ins| no-mount %.4f -> with mount %.4f N*m; ' ...
    'max |delta Torque_ins| over cells = %.4f N*m (deflection carried)\n'], ...
    m0, mA, dT);

% Material-ratio check (plan item 4 gate): stiffness scales with E and
% deflection with 1/E. Analytic on k (exact), approximate on the
% deflection (the solved force state re-solves slightly per material).
Eratio = 205 / 68.9;
kA_alu = 1 / oA.CmountEffective(1, 1);
cfgC = cfgA;
cfgC.mountMaterial = 'chromoly4130';
oC = build(cfgC);
kA_cro = 1 / oC.CmountEffective(1, 1);
kRatio = kA_cro / kA_alu;
dC = maxDefl(oC.MountDeflection(:, :, :, aK));
dRatio = dA / dC;
fprintf('  k(axial) aluminum = %.4g N/m, CroMoly = %.4g N/m; ratio = %.9f (E ratio %.9f)\n', ...
    kA_alu, kA_cro, kRatio, Eratio);
fprintf('  max deflection aluminum/CroMoly = %.4f (expect ~ %.4f)\n', dRatio, Eratio);
assert(abs(kRatio - Eratio) / Eratio < 1e-12, 'A: stiffness does not scale with E')
assert(abs(dRatio - Eratio) / Eratio < 0.10, ...
    'A: deflection ratio %.3f off the E ratio %.3f by >10%%', dRatio, Eratio)

% PENDING-BEN default-section behavior (never silent): material without a
% section prints the stand-in banner and still resolves a compliance.
fprintf('  --- material-without-section banner demo ---\n');
cfgP = S.spec.pulley;
cfgP.mountMaterial = 'aluminum6061T6';    % no mountSection on purpose
oP = build(cfgP);
assert(any(oP.CmountEffective(:)), 'A: PENDING-BEN default section gave zero compliance')

%% ---------------- Build B: bowden anchor deflection ------------------
fprintf('=== Build B: bowden anchor deflection (CroMoly 4130) ===\n');
cfgB = S.spec.pulley;
cfgB.routingMode = 'bowden';
cfgB.bowdenBracketMaterial = 'chromoly4130';
cfgB.bowdenBracketSection = secStandIn;
oB = build(cfgB);
[N1, N2] = oB.gridSize();

for k = 1:2
    dBk = oB.BowdenDeflection(:, :, :, k);    %N1 x 3 x N2 slice
    assert(any(~isnan(dBk(:))), 'B: BowdenDeflection all-NaN at crossing %d', k)
    fin = dBk(~isnan(dBk));
    assert(all(isfinite(fin)), 'B: BowdenDeflection not finite at crossing %d', k)
    % u_t pinned at the theta0 cell in the child frame (base bowden rule).
    ii0 = mod(oB.RefIndex(k) - 1, N1) + 1;
    iii0 = floor((oB.RefIndex(k) - 1) / N1) + 1;
    vref = reshape(oB.u_t(ii0, :, iii0, k), 1, 3, 1);   %components on dim 2
    dev = max(abs(oB.u_t(:, :, :, k) - repmat(vref, [N1, 1, N2])), [], 'all');
    assert(dev == 0, 'B: bowden u_t not pinned at crossing %d (dev %g)', k, dev)
    fprintf('  crossing %d: pinned u_t (max dev %.1e), max |anchor deflection| = %.4f mm\n', ...
        k, dev, 1000 * maxDefl(dBk));
end
mB = minTau(oB, oB.ActiveCrossing);
fprintf('  min |Torque_ins| (bowden, active crossing %d) = %.4f N*m\n', ...
    oB.ActiveCrossing, mB);

% Packaging envelope: Shimano-type defaults (scalar when the config omits
% the fields -- cfgField2 broadcasts only supplied values) + tunable.
assert(isequal(oB.BowdenBossDia, 0.007) && ...
    isequal(oB.BowdenRunClearance, 0.005), ...
    'B: default bowden envelope fields wrong')
cfgT = cfgB;
cfgT.bowdenBossDia = 0.008;         %tunability check
cfgT.bowdenRunClearance = 0.006;
oT = build(cfgT);
assert(isequal(oT.BowdenBossDia, [0.008, 0.008]) && ...
    isequal(oT.BowdenRunClearance, [0.006, 0.006]), ...
    'B: bowden envelope fields not tunable')
fprintf('  envelope defaults boss 7.0 mm / clearance 5.0 mm; tuned 8.0 / 6.0 mm OK\n');

fprintf('=== live_builds_BiPulley2_20261022: ALL CHECKS PASS ===\n');

function d = maxDefl(grid3)
% max cell magnitude of an N1x3xN2(x2) deflection grid over the 3-axis;
% NaN cells ("no bracket here" slots) count as zero.
n = vecnorm(grid3, 2, 2);
n(isnan(n)) = 0;
d = max(n(:));
end

function m = minTau(o, k)
tauAbs = abs(o.Torque_ins(:, :, :, k));
validTau = tauAbs(~isnan(tauAbs));
if isempty(validTau)
    m = NaN;
else
    m = min(validTau);
end
end

function d = maxDeltaTau(oa, ob, k)
% Max |Torque_ins| change between two builds over jointly non-NaN cells.
a = oa.Torque_ins(:, :, :, k);
b = ob.Torque_ins(:, :, :, k);
m = ~isnan(a) & ~isnan(b);
d = max(abs(a(m) - b(m)));
end
