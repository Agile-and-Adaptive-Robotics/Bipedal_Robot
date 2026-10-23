% Verify (2026-10-22): deterministic re-grade of the smoke's saved best
% design. Loads the rung-2 mat, rebuilds the graded grid from its own
% spec/T, re-evaluates xBest through BiPam_pulley2_20261022, recomputes
% the two percentages against the human CSV, and asserts they equal the
% recorded surveyMeta values. Deterministic (class construction +
% arithmetic, no RNG).
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));

S = load(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_survey_tib_ant_r_nB2_G1_20261022_133409.mat'));
m = S.surveyMeta;
fprintf('mat meta: %s | nB%d G%d | meet %.2f%% active %.2f%%\n', ...
    m.seed, m.nBPA, m.G, m.meetPct, m.activePct);

obj = BiPam_pulley2_20261022(S.spec.name, S.points, S.spec.cross, 20, ...
    S.T, S.xBest(1), 0.75 * S.xBest(1), S.xBest(2), 0.025, 620, ...
    0.005, 5e5, 1e4, 6, S.spec.pulley);
axisU = S.spec.kin.axes(m.slot, :) / norm(S.spec.kin.axes(m.slot, :));
[N1, N2] = obj.gridSize();
T3 = obj.Torque_ins(:, :, :, m.slot);
tauB = reshape(sum(T3 .* reshape(axisU, 1, 3, 1), 2), N1, N2);
strain = reshape(obj.strain_p, N1, N2);
slack = logical(obj.PulleySlack);
infeas = logical(obj.PulleyInfeasible);
taut = ~slack & ~infeas & isfinite(tauB);
tauH = S.tauHuman;
meet = taut & isfinite(tauH) & ...
    (sign(tauB) == sign(tauH) | abs(tauH) <= 1e-9) & ...
    abs(tauB) >= abs(tauH) - 1e-12;
active = taut & isfinite(strain) & strain > 0 & abs(tauB) > 0;
meetPct = 100 * sum(meet(:)) / numel(meet);
activePct = 100 * sum(active(:)) / numel(active);
fprintf('re-grade: meet %.2f%% active %.2f%% (grid %dx%d)\n', ...
    meetPct, activePct, N1, N2);
assert(abs(meetPct - m.meetPct) < 1e-9, 'meet mismatch');
assert(abs(activePct - m.activePct) < 1e-9, 'active mismatch');
fprintf('VERIFY PASS: recorded percentages reproduced exactly.\n');
