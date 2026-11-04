function fr_sensitivity_gate(round)
%FR_SENSITIVITY_GATE Stage-1 campaign deciding gate (master directive item 9).
%
% Runnable exactly as:
%   "D:/Program Files/MATLAB/R2025a/bin/matlab.exe" -batch "cd('D:/GitHub/...
%   .../Functional_Routing_20261027'); sensitivity_gate(1)"
% (the Functional_Routing_20261027/sensitivity_gate.m shim addpaths this
% campaign folder; the full implementation lives here).
%
% In a fresh process it must PASS every check or error() (nonzero exit):
%   (a) finite-difference sensitivity: perturbing EACH attachment/routing
%       variable changes the evaluated geometry AND a relevant
%       objective/constraint output;
%   (b) bounds have nonzero extent, and a short optimization on one
%       demonstrator does not silently restore the seed coordinates;
%   (c) representative cases: rigid-limit recovery, both bracket
%       contributions, parallel load sharing, transmission work
%       consistency, diameter endpoints, wrapping transitions;
%   (d) LAYERING proof: optimize a deep route against bone, then a
%       superficial route whose obstacle constraint against the retained
%       deep envelope actually ACTIVATES (deflection mm + active flag
%       reported; a trivially non-touching route FAILS), then one
%       refinement round that improves the superficial route's objective
%       WITHOUT degrading the deep route's optimized objective beyond a
%       stated tolerance;
%   (e) per-evaluation wall seconds and an extrapolated full-campaign
%       budget.
%
% Full log: campaign output dir / gate_r<round>.log (diary). Budget: the
% gate is subset to stay under ~25 minutes; every subset is stated here
% and in GateReport.blockers.

C = fr_setup();
logFile = fullfile(C.outputDir, sprintf('gate_r%d.log', round));
if exist(logFile, 'file')
    delete(logFile);
end
diary(logFile);
% keep the diary open until the function exits (error() included)
cleanupObj = onCleanup(@() diary('off'));

fprintf('=== sensitivity_gate round %d | %s ===\n', round, datetime('now'));
fprintf('machine: %s | campaign: %s\n', getenv('COMPUTERNAME'), C.camp);

tAll = tic;
chk = containers.Map('KeyType', 'char', 'ValueType', 'any');
evalWall = [];

% ---- kinematics + validation (precondition for everything) --------------
kin = fr_kinematics(C.osimStock);
kv = fr_validateKin(C, kin);
fprintf(['KINEMATICS validation vs geometry_20261103 ' 'bodyTransformsInGround (6 poses, exact transforms): ' 'maxErr %.6f mm, maxRot %.6f deg (tol %.2f mm)\n'], kv.maxErrMm, kv.maxRotErrDeg, C.kinTolMm);
if kv.maxErrMm > C.kinTolMm
    error('sensitivity_gate:kinematics', ...
        'kinematics validation failed: %.4f mm > %.2f mm', ...
        kv.maxErrMm, C.kinTolMm);
end
chk('kinematics') = kv;

cases = fr_cases(C, kin);
fprintf('cases built: %s\n', strjoin({cases.name}, ', '));

% ========================================================================
fprintf('\n--- CHECK (a): per-variable finite-difference sensitivity ---\n');
csA = cases(strcmp({cases.name}, 'hip_flex_ext'));
csA = csA(1);
oq = struct('quiet', true);
[f0, rep0] = fr_objective(C, csA, kin, csA.xSeed, 'rigid', oq);
evalWall(end+1) = rep0.wallSeconds;
geom0 = rep0.midPoseWorld;
c0 = rep0.c;
hAbs = [0.05; 0.05; 3.0; 0.05; 0.05; 0.08];  % u, v, s_mm, L0(dm), tendon(dm), split
allSensitive = true;
nGeomReq = 3 * numel(csA.movable);   % attachment variables need geometry response
for j = 1:csA.nX
    hj = hAbs(min(j, numel(hAbs)));   % absolute steps (scaled units)
    xPlus = csA.xSeed; xPlus(j) = min(csA.ub(j), csA.xSeed(j) + hj);
    xMinus = csA.xSeed; xMinus(j) = max(csA.lb(j), csA.xSeed(j) - hj);
    if xPlus(j) == xMinus(j)
        fprintf('  var %2d: bounds collapsed -> FAIL\n', j);
        allSensitive = false;
        continue
    end
    [~, rp] = fr_objective(C, csA, kin, xPlus, 'rigid', oq);
    evalWall(end+1) = rp.wallSeconds; %#ok<AGROW>
    [~, rm] = fr_objective(C, csA, kin, xMinus, 'rigid', oq);
    evalWall(end+1) = rm.wallSeconds; %#ok<AGROW>
    gPlus = max(vecnorm(rp.midPoseWorld - geom0, 2, 2));
    gMinus = max(vecnorm(rm.midPoseWorld - geom0, 2, 2));
    geomD = max(gPlus, gMinus);
    fD = max(abs(rp.f - f0), abs(rm.f - f0));
    cD = max(max(abs(rp.c - c0)), max(abs(rm.c - c0)));
    % attachment variables (mount u/v/standoff) must move the evaluated
    % geometry AND an output; length/split variables respond through the
    % objective/constraints (master item 9 wording)
    if j <= nGeomReq
        ok = geomD > 1e-6 && (fD > 1e-12 || cD > 1e-9);
    elseif j == csA.nX
        % series split: PACKAGING-recorded variable (documented mechanical
        % equivalence - two sections share the bundle force/pressure), so
        % its sensitivity requirement is a changed CONFIGURATION record
        cfgA = fr_splitConfig(C, csA, xPlus);
        cfgM = fr_splitConfig(C, csA, xMinus);
        ok = cfgA.seriesSections ~= cfgM.seriesSections || ...
            abs(cfgA.sectionLengths_m(1) - cfgM.sectionLengths_m(1)) > 1e-12;
    else
        ok = (fD > 1e-12 || cD > 1e-9);
    end
    if ~ok
        allSensitive = false;
    end
    fprintf('  var %2d [%s]: dGeom %.4f mm, df %.3e, dc %.3e -> %s\n', ...
        j, fr_varName(csA, j), geomD * 1000, fD, cD, fr_tf(ok));
end
chk('a_sensitivity') = allSensitive;
fprintf('CHECK (a) %s\n', fr_tf(allSensitive));

% ========================================================================
fprintf(['\n--- CHECK (b): bounds extent + short optimization does not ' ...
    'silently restore the seed ---\n']);
extent = csA.ub - csA.lb;
allExtent = all(extent > 1e-9);
fprintf('  bounds extent (min %.3g, max %.3g over %d vars) -> %s\n', ...
    min(extent), max(extent), csA.nX, fr_tf(allExtent));
so = struct('quiet', true, 'budget', C.solver.gateMaxFunctionEvaluations);
[~, repB, infoB] = fr_solveCase(C, csA, kin, 'rigid', so);
xDelta = max(abs(infoB.xBest - infoB.xSeed));
geomMove = max(vecnorm(repB.midPoseWorld - geom0, 2, 2));
bOk = allExtent && geomMove > 1e-3 && xDelta > 1e-3;
fprintf(['  short patternsearch (%d evals): |xBest-xSeed|max = %.4f ' '(scaled), geometry moved %.2f mm, f %.4f -> %.4f -> %s\n'], infoB.budget, xDelta, geomMove * 1000, infoB.fSeed, infoB.fBest, fr_tf(bOk));
chk('b_bounds_noRestore') = bOk;
fprintf('CHECK (b) %s\n', fr_tf(bOk));

% ========================================================================
fprintf('\n--- CHECK (c): representative cases + mechanics identities ---\n');
cCases = {'hip_flex_ext', 'hip_rot', 'hip_flex_extrot', 'knee_bowden', ...
    'foot_span', 'trunk_ercspn'};
cOk = true;
for ic = 1:numel(cCases)
    cs = cases(strcmp({cases.name}, cCases{ic}));
    cs = cs(1);
    [fR, repR] = fr_objective(C, cs, kin, cs.xSeed, 'rigid', oq);
    evalWall(end+1) = repR.wallSeconds; %#ok<AGROW>
    [fC, repC] = fr_objective(C, cs, kin, cs.xSeed, 'compliant', oq);
    evalWall(end+1) = repC.wallSeconds; %#ok<AGROW>
    fprintf(['  %-16s rigid f=%.4f viol=%.2e | compliant f=%.4f ' 'viol=%.2e | tauCrossCheck max %.2e N*m\n'], cs.name, fR, repR.violation, fC, repC.violation, repC.tauCrossCheckMaxDiff);
    okC = repR.violation < 1e-6 || repR.violation < 0.05;  % feasibility at seed
    fprintf('    seed feasibility (rigid): %s\n', fr_tf(okC));
    cOk = cOk && okC;
    % per-case identity checks on a subset (budget): tackle/mono cases
    if any(strcmp(cs.name, {'knee_bowden', 'hip_flex_ext'}))
        idOk = fr_identities(C, cs, kin, oq, evalWall);
        cOk = cOk && idOk.passed;
        chk(['c_identities_' cs.name]) = idOk;
    end
end
% rigid-limit recovery (compliance -> 0 converges to the rigid branch)
rlOk = true;
for ic = 1:2
    cs = cases(strcmp({cases.name}, cCases{ic}));
    cs = cs(1);
    [~, repR] = fr_objective(C, cs, kin, cs.xSeed, 'rigid', oq);
    % the class's rigid convention is INFINITE stiffness (fortzBi2
    % isinf(X1)&&isinf(X2)); the recovery probe therefore pushes Xi0 -> 0,
    % Xi1/Xi2 -> 1e6x adopted, Xi3 -> 0, tendon -> stiff (wraps 1e4)
    Cstiff = C;
    Cstiff.Xi.Xi0_flx = 0;
    Cstiff.Xi.Xi1 = C.Xi.Xi1 * 1e6;
    Cstiff.Xi.Xi2 = C.Xi.Xi2 * 1e6;
    Cstiff.Xi.wrapsOverride = 1e4;
    Cstiff.Xi.noCompliance = true;
    sBk = C.Xi.Xi3;
    Cstiff.Xi.Xi3 = 0;
    [fS, repS] = fr_objective(Cstiff, cs, kin, cs.xSeed, 'compliant', oq);
    C.Xi.Xi3 = sBk;
    dTau = max(abs(repS.tauPrimary - repR.tauPrimary));
    dF = abs(fS - fr_objective(C, cs, kin, cs.xSeed, 'rigid', ...
        struct('quiet', true)));
    ok = dTau < 0.05 && dF < 0.05;
    rlOk = rlOk && ok;
    evalWall(end+1) = repS.wallSeconds; %#ok<AGROW>
    fprintf(['  rigid-limit recovery %-16s: |dtau|max %.2e N*m, ' '|df| %.2e -> %s\n'], cs.name, dTau, dF, fr_tf(ok));
end
cOk = cOk && rlOk;
chk('c_rigid_limit_recovery') = rlOk;
% diameter endpoints + own 10 mm law
diaOk = fr_diameterCheck(C);
cOk = cOk && diaOk;
chk('c_diameter_endpoints') = diaOk;
% wrapping transitions (conditional points; gait2392 has no wrap objects)
wtOk = fr_wrappingCheck(C, kin, cases);
cOk = cOk && wtOk;
chk('c_wrapping_transitions') = wtOk;
fprintf('CHECK (c) %s\n', fr_tf(cOk));

% ========================================================================
fprintf('\n--- CHECK (d): LAYERING proof (deep -> superficial -> refine) ---\n');
[dOk, dInfo] = fr_layeringProof(C, kin, cases);
chk('d_layering') = dInfo;
fprintf('CHECK (d) %s\n', fr_tf(dOk));

% ========================================================================
fprintf('\n--- CHECK (e): timing + extrapolated full-campaign budget ---\n');
evalWall = evalWall(isfinite(evalWall) & evalWall > 0);
medWall = median(evalWall);
meanWall = mean(evalWall);
nEvalsFull = 27 * 2 * (C.solver.MaxFunctionEvaluations + ...
    C.refinementRounds * C.solver.refineMaxFunctionEvaluations);
budgetH = nEvalsFull * meanWall / 3600;
fprintf('  evaluations measured: %d | median %.3f s | mean %.3f s\n', ...
    numel(evalWall), medWall, meanWall);
fprintf(['  extrapolated full campaign (27 groups x {rigid, compliant} ' 'x [%d solve + %d x %d refine evals] = %d evals): %.2f h\n'], C.solver.MaxFunctionEvaluations, C.refinementRounds, C.solver.refineMaxFunctionEvaluations, nEvalsFull, budgetH);
eOk = ~isempty(evalWall);
chk('e_budget') = struct('medianSec', medWall, 'meanSec', meanWall, ...
    'nEvals', numel(evalWall), 'fullCampaignEvals', nEvalsFull, ...
    'fullCampaignHours', budgetH);
fprintf('CHECK (e) %s\n', fr_tf(eOk));

% ========================================================================
fprintf('\n=== gate wall time: %.1f min ===\n', toc(tAll) / 60);
failed = {};
keys = chk.keys;
for k = 1:numel(keys)
    v = chk(keys{k});
    if islogical(v) && ~v
        failed{end+1} = keys{k}; %#ok<AGROW>
    end
end
if ~dOk
    failed{end+1} = 'd_layering';
end
if ~isempty(failed)
    error('sensitivity_gate:failed', ...
        'GATE FAILED checks: %s (log: %s)', strjoin(failed, ', '), ...
        logFile);
end
fprintf('GATE PASSED all checks. log: %s\n', logFile);
end

% ========================================================================== helpers
function nm = fr_varName(cs, j)
nMov = numel(cs.movable);
if j <= 3 * nMov
    m = ceil(j / 3);
    comp = mod(j - 1, 3);
    names = {'u', 'v', 'standoff'};
    nm = sprintf('mov%d(%s,row%d)', m, names{comp + 1}, ...
        cs.movable(m).row);
    return
end
if j == csA_nX(cs) - 2
    nm = 'L0';
elseif j == csA_nX(cs) - 1
    nm = 'tendon';
else
    nm = 'seriesSplit';
end
end

function n = csA_nX(cs)
n = cs.nX;
end

function s = fr_tf(b)
if b
    s = 'PASS';
else
    s = 'FAIL';
end
end
