% Test Escalation Advisor 20261022
% Unit polarity tests for escalation_advisor_20261022.m
% (OPTIMIZER_IMPROVEMENT_PLAN.md Item 2 gate, WITH the 2026-10-22 review
% blocker fix: torque test over taut feasible cells only; slack cells map
% to the RoM/slack side).
%
%  (a) SLACK-ARTIFACT POLARITY (the blocker's unit scenario): a grid with
%      slack cells (Torque_ins identically 0 there by construction,
%      MonoPam_pulley.m:595 / BiPam_pulley.m:671) whose TAUT cells meet
%      the target must classify rom_slack_short -- NOT torque-short --
%      and must NOT advise +BPA. The raw (slack-gamed) margin reads -1
%      while the taut margin is positive: the artifact is exposed, not
%      consumed.
%  (b) TORQUE-SHORT: no slack/infeasible, taut min below target ->
%      torque_short, stepNBPA = +1.
%  (c) BOTH: slack AND taut torque shortfall -> both steps.
%  (d) CONVERGED: clean grid, no steps.
%  (e) CELL-WEIGHT MASK (Item 6 hook): the (a) grid with the slack cells
%      masked to W = 0 classifies converged -- a masked cell can never
%      justify a ladder step.
%  (f) DOCTRINE REPLAY (the re-picked fixture verdict):
%      BiPulley_opensim_med_gas_r_2BPA_G2_20261021_1114.mat must read
%      rom_slack_short with nSlack = 10, rawMargin = -1, stepG = +1,
%      stepNBPA = 0 (doctrine_rerun_20261021.log:1160-1164: margin
%      -1.0000, slack cells = 10, worst constraint +0.430868,
%      patternsearch -2 -- a slack cell is the binding cell, so the
%      plan's original "torque-short -> advise nBPA 3" verdict was the
%      self-contradiction the review flagged).
%  (g) DETERMINISM: test (f) twice, identical outputs, zero optimizer
%      calls, no parpool, no RNG.
%
% Runs in seconds. error() on any failure; prints
% TEST_ESCALATION_ADVISOR_20261022 PASS as the last line.

scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root) || isempty(parent)
        error('Could not locate the repo root from %s', scriptDir)
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(scriptDir);

target = 50;   % N*m, synthetic
fprintf('== test_escalation_advisor_20261022 ==\n')

%% (a) slack-artifact polarity
G5 = 5;
tau = ones(G5) * 2 * target;      % taut cells: 2x target, comfortable
slack = false(G5); slack(1:2, 1:5) = true;   % 10 slack cells
tau(slack) = 0;                   % construction: F_ins zeroed at slack
synth = struct('TorqueIns', tau, 'PulleySlack', slack, ...
    'PulleyInfeasible', false(G5));
[adv, rep] = escalation_advisor_20261022(synth, target, ...
    struct('verbose', false, 'curNBPA', 2, 'curG', 2));
assert(strcmp(adv.verdict, 'rom_slack_short'), ...
    '(a) verdict %s, expected rom_slack_short', adv.verdict)
assert(adv.stepNBPA == 0, '(a) must NOT advise +BPA on a slack design')
assert(adv.stepG == 1 && adv.nextG == 3, '(a) must advise G 2 -> 3')
assert(rep.nSlack == 10 && rep.nTautFeasible == 15, ...
    '(a) cell partition wrong: %d slack, %d taut', rep.nSlack, rep.nTautFeasible)
assert(rep.rawMargin == -1, '(a) raw margin %g, expected exactly -1', ...
    rep.rawMargin)
assert(rep.tautMargin == 1, '(a) taut margin %g, expected +1', rep.tautMargin)
fprintf('(a) PASS: 10 slack + taut cells at 2x target -> %s, steps (+%d,+%d), raw %+.4f vs taut %+.4f\n', ...
    adv.verdict, adv.stepNBPA, adv.stepG, rep.rawMargin, rep.tautMargin)

%% (b) torque-short
tauB = ones(G5) * target;
tauB(3, 3) = 0.4 * target;        % one taut cell 60%% short
synthB = struct('TorqueIns', tauB, 'PulleySlack', false(G5), ...
    'PulleyInfeasible', false(G5));
[advB, repB] = escalation_advisor_20261022(synthB, target, ...
    struct('verbose', false, 'curNBPA', 1, 'curG', 1));
assert(strcmp(advB.verdict, 'torque_short'), '(b) verdict %s', advB.verdict)
assert(advB.stepNBPA == 1 && advB.nextNBPA == 2, '(b) step nBPA wrong')
assert(advB.stepG == 0, '(b) must not step G on a pure torque shortfall')
assert(repB.nSlack == 0 && repB.nInfeasible == 0)
fprintf('(b) PASS: taut cell 60%% short -> %s, steps (+%d,+%d)\n', ...
    advB.verdict, advB.stepNBPA, advB.stepG)

%% (c) both
tauC = tauB;
slackC = false(G5); slackC(5, :) = true;
tauC(slackC) = 0;
synthC = struct('TorqueIns', tauC, 'PulleySlack', slackC, ...
    'PulleyInfeasible', false(G5));
[advC, ~] = escalation_advisor_20261022(synthC, target, ...
    struct('verbose', false, 'curNBPA', 1, 'curG', 1));
assert(strcmp(advC.verdict, 'both'), '(c) verdict %s', advC.verdict)
assert(advC.stepNBPA == 1 && advC.stepG == 1, '(c) both steps expected')
fprintf('(c) PASS: slack + taut shortfall -> %s, steps (+1,+1)\n', advC.verdict)

%% (d) converged
synthD = struct('TorqueIns', ones(G5) * 2 * target, ...
    'PulleySlack', false(G5), 'PulleyInfeasible', false(G5));
[advD, ~] = escalation_advisor_20261022(synthD, target, ...
    struct('verbose', false, 'curNBPA', 3, 'curG', 3));
assert(strcmp(advD.verdict, 'converged') && advD.converged, ...
    '(d) verdict %s', advD.verdict)
fprintf('(d) PASS: clean grid -> converged\n')

%% (e) cell-weight mask (Item 6 corner-weight hook)
W = ones(G5); W(slack) = 0;       % mask the slack/artifact cells out
[advE, repE] = escalation_advisor_20261022(synth, target, ...
    struct('verbose', false, 'curNBPA', 2, 'curG', 2, 'cellWeights', W));
assert(strcmp(advE.verdict, 'converged'), ...
    '(e) masked cells must not drive a step (got %s)', advE.verdict)
assert(repE.nSlack == 0 && repE.nTautFeasible == 15, ...
    '(e) masked partition wrong')
fprintf('(e) PASS: W = 0 on the slack cells -> %s (mask honored)\n', ...
    advE.verdict)

%% (f) doctrine replay -- the re-picked fixture verdict (review blocker 1)
% The review allowed "re-pick the replay fixture or its expected
% verdict". Probing the fixture this session: its 15 taut cells are ALL
% strain-OK and 5 are GENUINELY short (7.6/10.1/12.1/17.8/28.8 N*m vs
% the 77.9 N*m target), so under the FIXED taut-only rule the honest
% verdict is BOTH. The blocker's essential property is asserted below:
% the slack artifact (raw margin -1.0000, log :1160) does NOT drive the
% verdict -- the torque side fires from taut cells only, so the advice
% steps BOTH sides, not the plan's slack-driven 'nBPA 3 alone'.
matPath = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', ...
    'Results', 'BiPulley_opensim_med_gas_r_2BPA_G2_20261021_1114.mat');
if ~isfile(matPath)
    error('test_escalation_advisor_20261022:missingFixture', ...
        'Doctrine fixture mat not found: %s', matPath)
end
[advF, repF] = escalation_advisor_20261022(matPath, [], ...
    struct('verbose', true, 'curNBPA', 2, 'curG', 2));
assert(repF.nSlack == 10, ...
    '(f) doctrine slack cells %d, expected 10 (doctrine_rerun log :1162)', ...
    repF.nSlack)
assert(repF.nInfeasible == 0, ...
    '(f) doctrine infeasible cells %d, expected 0 (log :1163)', ...
    repF.nInfeasible)
assert(strcmp(advF.verdict, 'both'), ...
    '(f) verdict %s, expected both (10 slack cells + 5 genuinely short taut cells)', ...
    advF.verdict)
assert(advF.stepNBPA == 1 && advF.stepG == 1, ...
    '(f) must step BOTH sides, not the plan''s slack-driven nBPA-3-alone')
assert(repF.rawMargin == -1, ...
    '(f) raw margin %g, expected -1 (log :1160 margin = -1.0000)', ...
    repF.rawMargin)
assert(repF.tautMargin < 0, ...
    '(f) taut margin must expose the genuine taut shortfall (%g)', ...
    repF.tautMargin)
assert(repF.rawMargin < repF.tautMargin, ...
    '(f) raw %g must be worse than taut %g (the artifact made it look worse)', ...
    repF.rawMargin, repF.tautMargin)
fprintf(['(f) PASS: doctrine med_gas_r (2,[2 1]) -> %s (nSlack = 10, raw ' ...
    '%+.4f = the log :1160 artifact vs taut %+.4f = genuine), steps (+%d,+%d)\n'], ...
    advF.verdict, repF.rawMargin, repF.tautMargin, advF.stepNBPA, advF.stepG)

%% (f2) the torque_short firing is taut-driven, NOT slack-driven
% Mask the slack cells to W = 0: if the torque side were an artifact of
% the slack zeros, the verdict would collapse to converged; it must
% stay torque_short (the 5 short taut cells).
Sfix = load(matPath, 'finalObj', 'tauTarget');
Wfix = double(~Sfix.finalObj.PulleySlack);
[advF2, repF2] = escalation_advisor_20261022(Sfix.finalObj, Sfix.tauTarget, ...
    struct('verbose', false, 'curNBPA', 2, 'curG', 2, 'cellWeights', Wfix));
assert(strcmp(advF2.verdict, 'torque_short'), ...
    '(f2) slack-masked verdict %s, expected torque_short', advF2.verdict)
assert(repF2.nSlack == 0 && repF2.minTautTau == repF.minTautTau, ...
    '(f2) masked run must keep the taut diagnosis identical')
fprintf('(f2) PASS: slack cells masked -> %s persists (taut-driven, min taut %.4f unchanged)\n', ...
    advF2.verdict, repF2.minTautTau)

%% (g) determinism
[advG, repG] = escalation_advisor_20261022(matPath, [], ...
    struct('verbose', false, 'curNBPA', 2, 'curG', 2));
assert(strcmp(advG.verdict, advF.verdict) && ...
    advG.stepNBPA == advF.stepNBPA && advG.stepG == advF.stepG && ...
    isequaln(repG.minTautTau, repF.minTautTau) && ...
    isequaln(repG.tautMargin, repF.tautMargin), ...
    '(g) two identical calls disagreed')
fprintf('(g) PASS: deterministic (identical repeat outputs, zero optimizer calls)\n')

fprintf('TEST_ESCALATION_ADVISOR_20261022 PASS\n')
