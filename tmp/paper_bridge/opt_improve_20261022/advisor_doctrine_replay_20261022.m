% Advisor Doctrine Replay 20261022
% Replays escalation_advisor_20261022 over the 2026-09-26 G = 1 baseline
% mats AND the 2026-10-21 doctrine mats (read-only; no optimizer calls,
% no parpool). Evidence run for the review blocker-1 fix: the plan's
% named fixture (med_gas_r at (2,[2 1]), doctrine_rerun_20261021.log:1160,
% margin -1.0000 with 10 slack cells) must classify rom_slack_short, NOT
% 'torque-short -> advise nBPA 3'.
% Writes tmp/paper_bridge/opt_improve_20261022/advisor_doctrine_replay_20261022.csv
% and prints one block per mat.

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

resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
d = [dir(fullfile(resDir, 'BiPulley_opensim_*_20260926_*.mat')); ...
    dir(fullfile(resDir, 'BiPulley_opensim_*_2BPA_G2_20261021_*.mat')); ...
    dir(fullfile(resDir, 'BiPulley_opensim_*_2BPA_G4_20261021_*.mat')); ...
    dir(fullfile(resDir, 'BiPulley_opensim_*_4BPA_G4_20261021_*.mat'))];
files = cellfun(@(f) fullfile(resDir, f), {d.name}, 'UniformOutput', false);
files = files(:);
fprintf('== advisor_doctrine_replay_20261022: %d anchor mats ==\n', numel(files))

csvPath = fullfile(scriptDir, 'advisor_doctrine_replay_20261022.csv');
if isfile(csvPath)
    delete(csvPath)   % this script OWNS its csv (not a record of record)
end
counts = struct('converged', 0, 'torque_short', 0, 'rom_slack_short', 0, 'both', 0);
for a = 1:numel(files)
    mNB = 1; mG = 1;   % ladder position defaults; read the mat's config
    try
        S = load(files{a}, 'spec');
        mNB = S.spec.pulley.nPulleyBPA;
        if isscalar(mNB), mNB = mNB; end
        tlp = S.spec.pulley.tackleLineParts;
        mG = max(tlp(:));
    catch
    end
    [advice, rep] = escalation_advisor_20261022(files{a}, [], ...
        struct('verbose', true, 'curNBPA', mNB, 'curG', mG, ...
        'csvFile', csvPath));
    if ismember(advice.verdict, fieldnames(counts))
        counts.(advice.verdict) = counts.(advice.verdict) + 1;
    end
    fprintf('  [%d/%d] %s -> verdict %-15s steps (+%d,+%d) | slack %d, infeas %d, taut %d | raw %+.4f taut %s\n', ...
        a, numel(files), rep.source, advice.verdict, advice.stepNBPA, ...
        advice.stepG, rep.nSlack, rep.nInfeasible, rep.nTautFeasible, ...
        rep.rawMargin, sprintf('%.4f', rep.tautMargin));
end
fprintf(['VERDICT COUNTS: converged %d | torque_short %d | rom_slack_short ' ...
    '%d | both %d\n'], counts.converged, counts.torque_short, ...
    counts.rom_slack_short, counts.both)
fprintf('CSV: %s\n', csvPath)
fprintf('ADVISOR_DOCTRINE_REPLAY_20261022 DONE\n')
