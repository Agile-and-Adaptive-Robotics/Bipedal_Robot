% replay_gate_BiPulley2_20261022 -- deterministic replay regression gate
% for the dated variant class BiPam_pulley2_20261022 (copied from
% BiPam_pulley.m; copy-to-new-variant, Ben's class untouched).
%
% GATE (the established 0.000-tolerance rebuild gate, doctrine round
% 2026-10-21 form; assembled here per the 2026-10-22 review blocker 2):
% at nBPA = 1, G = [1 1], routingMode 'moving_via' (the 2026-09-26
% recorded batch, Results/BiPulley_opensim_*_20260926_*.mat), re-evaluate
% each mat's recorded xBest through the VARIANT class with its DEFAULT
% config (mount/bowden compliances off) and assert the recomputed
% torques equal the recorded finalObj torques to 0.000 N*m.
% Deterministic by construction: evaluate() is festo4 + fzero, no RNG
% (Opt_run_BiPulley.m:360-385). The fixture's config is SELF-CHECKED
% before replay (lesson from review blocker 1: never trust a fixture's
% claimed configuration). Read-only over Results/.
%
% Run: matlab -batch "run('tmp/paper_bridge/opt_improve_20261022/replay_gate_BiPulley2_20261022.m')"

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

resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
mats = dir(fullfile(resDir, 'BiPulley_opensim_*_20260926_*.mat'));
assert(~isempty(mats), 'No 2026-09-26 batch mats found in %s', resDir);

TOL = 0.000;    % N*m -- the established gate
nFail = 0;
worstTau = -inf;
worstMargin = -inf;
worstState = -inf;

fprintf('=== replay_gate_BiPulley2_20261022: %d recorded mats ===\n', numel(mats));
for j = 1:numel(mats)
    f = fullfile(mats(j).folder, mats(j).name);
    S = load(f);

    % Fixture self-check: the gate is DEFINED at nBPA=1, G=[1 1],
    % moving_via (ask of 2026-10-22). A mat failing this is the wrong
    % fixture, not a gate failure.
    assert(isequal(S.finalObj.NPulleyBPA, [1, 1]), ...
        '%s: fixture nPulleyBPA is [%d %d], expected [1 1]', mats(j).name, ...
        S.finalObj.NPulleyBPA(1), S.finalObj.NPulleyBPA(2))
    assert(isequal(S.finalObj.PulleyGain, [1, 1]), ...
        '%s: fixture PulleyGain is [%g %g], expected [1 1]', mats(j).name, ...
        S.finalObj.PulleyGain(1), S.finalObj.PulleyGain(2))
    assert(all(strcmp(S.finalObj.RoutingMode, 'moving_via')), ...
        '%s: fixture RoutingMode is not moving_via/moving_via', mats(j).name)

    % evaluate() VERBATIM from Opt_run_BiPulley.m:360-385, class name
    % swapped to the variant. Default pulleyConfig -> variant compliances
    % off -> bit-identical code path.
    rest = S.xBest(1);
    tendon = S.xBest(2);
    kmax = 0.75 * rest;
    obj2 = BiPam_pulley2_20261022(S.spec.name, S.points, S.spec.cross, 20, ...
        S.T, rest, kmax, tendon, 0.025, S.basePressure, S.baseXi0, ...
        S.baseXi1, S.baseXi2, S.baseWraps, S.spec.pulley);
    aK = obj2.ActiveCrossing;
    tauAbs = abs(obj2.Torque_ins(:, :, :, aK));
    validTau = tauAbs(~isnan(tauAbs));
    margin2 = min(validTau) / S.tauTarget - 1;

    % Torque-grid compare (recorded batch finalObj torques), NaN pattern
    % must match too.
    A = obj2.Torque_ins(:);
    B = S.finalObj.Torque_ins(:);
    assert(isequal(isnan(A), isnan(B)), ...
        '%s: NaN pattern of Torque_ins differs from the record', mats(j).name)
    dTau = max(abs(A - B));
    dMargin = abs(margin2 - S.st.margin);

    % Whole-solve-state compare (stronger than the torque gate; the
    % recorded object carries the identical solve).
    dState = max([ ...
        max(abs(obj2.sContraction(:) - S.finalObj.sContraction(:))), ...
        max(abs(obj2.Ftendon(:) - S.finalObj.Ftendon(:))), ...
        max(abs(obj2.FsingleBPA(:) - S.finalObj.FsingleBPA(:))), ...
        max(abs(obj2.strain_p(:) - S.finalObj.strain_p(:))), ...
        max(abs(obj2.MuscleLength(:) - S.finalObj.MuscleLength(:))), ...
        1000 * double(~isequal(obj2.PulleySlack, S.finalObj.PulleySlack)), ...
        1000 * double(~isequal(obj2.PulleyInfeasible, S.finalObj.PulleyInfeasible))]);

    ok = (dTau <= TOL) && (dMargin <= TOL) && (dState <= 0);
    if ~ok
        nFail = nFail + 1;
    end
    fprintf(['%-46s dTorque_ins = %.6e N*m  dMargin = %.6e  ' ...
        'dSolveState = %.6e  [%s]\n'], mats(j).name, dTau, dMargin, ...
        dState, passFail(ok));
    worstTau = max(worstTau, dTau);
    worstMargin = max(worstMargin, dMargin);
    worstState = max(worstState, dState);
end
fprintf(['=== GATE %s: %d/%d mats reproduced | worst dTorque = %.6e N*m ' ...
    '(tol %.3f) | worst dMargin = %.6e | worst dSolveState = %.6e ===\n'], ...
    passFail(nFail == 0), numel(mats) - nFail, numel(mats), worstTau, ...
    TOL, worstMargin, worstState);
assert(nFail == 0, 'replay gate FAILED (%d mats)', nFail)

function s = passFail(ok)
if ok
    s = 'PASS';
else
    s = 'FAIL';
end
end
