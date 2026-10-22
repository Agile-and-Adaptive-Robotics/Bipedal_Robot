%TEST_SOFTCORNERWEIGHTS_20261022 Gate tests, plan Item 6 verification gate.
%   Deterministic, pure-function asserts; depends only on files in this
%   folder. Run from tmp/paper_bridge/opt_improve_20261022:
%     matlab -batch "test_softCornerWeights_20261022"
%   Gate items (plan lines 147-150):
%     * identity: with W = ones the weighted objective equals the
%       unweighted one bit-for-bit;
%     * polarity: a design failing ONLY in the tapered corner scores
%       strictly better than one failing mid-envelope; a mid-envelope
%       failure is never softened;
%     * record rule note: consumers report raw and softened margins.

fprintf('--- test_softCornerWeights_20261022 ---\n');
env1 = [-15 110];      % hip flexion task envelope (buildRoutingContext.m:75)
env2 = [-45 20];       % hip adduction graded box (actuator_map_27.json)
o = struct('band', 10, 'wMin', 0.1);   % PENDING-BEN defaults under test

% Identity: envelopes covering the whole grid => W == 1 exactly, and the
% weighted shortfall is bit-for-bit the unweighted one.
[t1, t2] = meshgrid(-20:2.5:120, -55:2.5:30);
Wide = softCornerWeights_20261022(t1, t2, [-1e3 1e3], [-1e3 1e3], o);
assert(isequal(Wide, ones(size(t1))), 'covering envelopes must give W == 1 exactly.');
S = 0.37 * mod(round(t1 * 3.7) + round(abs(t2) * 1.3), 7) / 7 - 0.2;  % deterministic
assert(isequal(Wide .* S, S), ...
    'with W == 1 the weighted shortfall must equal the raw one bit-for-bit.');
fprintf('PASS identity: W == 1 exactly; W .* S bit-for-bit == S\n');

% Mid-envelope never softened; floor exact at far corners
W = softCornerWeights_20261022(t1, t2, env1, env2, o);
assert(softCornerWeights_20261022(47.5, -12.5, env1, env2, o) == 1, ...
    'mid-envelope pose must carry W == 1 exactly.');
mid = t1 > 0 & t1 < 90 & t2 > -40 & t2 < 10;
assert(all(W(mid) == 1), 'every interior pose must carry W == 1 exactly.');
assert(softCornerWeights_20261022(130, -70, env1, env2, o) == o.wMin, ...
    'far corner must sit at the floor exactly.');
assert(softCornerWeights_20261022(-30, 40, env1, env2, o) == o.wMin, ...
    'opposite far corner must sit at the floor exactly.');
fprintf('PASS floor/mid: interior W == 1 exactly; far corners W == wMin == %g exactly\n', o.wMin);

% Polarity: corner-only failure scores strictly better than a mid failure
P = struct('band', 10, 'wMin', 0.1);
[t1p, t2p] = meshgrid(0:5:120, -40:5:25);
Wp = softCornerWeights_20261022(t1p, t2p, env1, env2, P);
[~, iA] = min(abs(t1p(:) - 105) + abs(t2p(:) - 25));   % tapered corner cell
[~, iB] = min(abs(t1p(:) - 45) + abs(t2p(:) - 0));    % mid-envelope cell
assert(Wp(iA) < 1 && Wp(iB) == 1, 'test cells are not in the intended weight zones.');
c = 0.5;                                              % shortfall magnitude
m1 = 0.05 * ones(size(Wp)); m1(iA) = -c;              % fails ONLY in the corner
m2 = 0.05 * ones(size(Wp)); m2(iB) = -c;              % fails ONLY mid-envelope
J1 = min(Wp(:) .* m1(:));                             % weighted worst margin
J2 = min(Wp(:) .* m2(:));
assert(J1 > J2, ...
    'corner-only failure (%g) must score strictly better than mid failure (%g).', J1, J2);
assert(Wp(iB) * m2(iB) == -c, 'the mid failure must be unchanged by weighting.');
fprintf(['PASS polarity: J(corner-only)=%g > J(mid-only)=%g; ' ...
    'mid failure untouched (== %g)\n'], J1, J2, -c);

% Composition with the Item 5 hip-coupling mask (AND via min): a pose that
% is INSIDE both joint envelopes but beyond the coupling limit gets caught.
aHip = [1 1];
[Whip] = hipAdductionEnvelope_20261022([92.5 92.5], [0 -45], ...
    struct('taperBandDeg', 10, 'wMin', 0.1));
Wtot = softCornerWeights_20261022([92.5 92.5], [0 -45], env1, env2, ...
    struct('band', 10, 'wMin', 0.1, 'hipMask', Whip));
assert(Whip(1) == 1 && Wtot(1) == 1, 'coupled-inside pose must stay fully trusted.');
assert(Whip(2) < 1, 'coupling mask must down-weight the combined-extreme pose.');
assert(Wtot(2) == min(1, Whip(2)) && Wtot(2) < 1, ...
    'hip mask must AND into the corner weights.');
fprintf('PASS composition: W(92.5,-45): corners-only=%g, hip-only=%.4f, AND=%.4f\n', ...
    softCornerWeights_20261022(92.5, -45, env1, env2, P), Whip(2), Wtot(2));

% Record rule (gate item 3) and slack guardrail reminder (review blocker):
fprintf(['NOTE record rule: consumers must report BOTH worstMarginRaw and ' ...
    'worstMarginSoft (never quote a softened number as raw).\n']);
fprintf(['NOTE slack guardrail: slack cells (F_ins zeroed: MonoPam_pulley.m:595, ' ...
    'BiPam_pulley.m:671) are EXCLUDED by mask from torque tests and routed to ' ...
    'the RoM/slack side; these weights never rescue them.\n']);

fprintf('TEST softCornerWeights_20261022: ALL PASS\n\n');
