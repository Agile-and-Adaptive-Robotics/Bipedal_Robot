function fr_perMuscleCompare()
%FR_PERMUSCLECOMPARE Robot BPA torque vs each NAMED human counterpart
%(Ben's 2026-11-04 rerun-mode request): instead of scoring each route
% only against its GROUP's summed signed target curve, this tool compares
% the optimized robot route's torque curve against EACH constituent
% human muscle's own signed per-pose curve (already stored per cell in
% the targets_v2 JSONs by build_signed_targets.py - no OpenSim rerun).
%
% For every group with run_1 checkpoints and every branch:
%   - rebuild the winner case (sizing_winners.json overrides)
%   - evaluate the design on the optimization grid
%   - per constituent muscle m: meetFrac_m = fraction of cells with
%     |tau_robot| >= |tau_m|; maxDeficit_m in N*m
%   - group-sum view (the scoring target of record) alongside
% Output: output/per_muscle_compare.csv (group, branch, muscle, ...) and
% output/per_muscle_compare.log. Rerun-mode OPTIMIZATION (making the
% route meet each named muscle individually) is a Ben-ruling-gated
% variant: it is NOT wired into the objective - this tool only measures.
%
% Run headless: matlab -batch "addpath('<camp>'); fr_perMuscleCompare"

logFile = fullfile('output', 'per_muscle_compare.log');
if ~exist('output', 'dir'); mkdir('output'); end
fid = fopen(logFile, 'w');
cleanupFid = onCleanup(@() fclose(fid));
fprintf(fid, 'per-muscle comparison | %s\n', datetime('now'));

C = fr_setup();
if isfield(C, 'kinematics'); kin = C.kinematics; else; kin = fr_kinematics(C.osimStock); end
szFile = fullfile(C.camp, 'sizing_winners.json');
assert(exist(szFile, 'file') == 2, 'fr_perMuscleCompare:noSizing', ...
    'sizing_winners.json missing - run fr_sizingSearch first');
sz = jsondecode(fileread(szFile));
if iscell(sz.groups); groupsList = sz.groups;
else; groupsList = arrayfun(@(g) g, sz.groups, 'UniformOutput', false);
end
ck1 = fullfile(C.camp, 'run_1', 'checkpoints');

csvFile = fullfile('output', 'per_muscle_compare.csv');
fcsv = fopen(csvFile, 'w');
fprintf(fcsv, ['group,branch,muscle,meetFrac,maxDeficit_Nm,' ...
    'peakAbsTauMuscle_Nm,peakAbsTauRobot_Nm\n']);

nEval = 0;
for gi = 1:numel(groupsList)
    gsz = groupsList{gi};
    if ~isfield(gsz, 'diaStr'); continue; end
    % winner case at the winner config
    ov = struct('diaStr', gsz.diaStr, 'nBPA', gsz.nBPA, ...
        'family', gsz.family, 'gain', gsz.gain);
    om = containers.Map('KeyType', 'char', 'ValueType', 'any');
    om(gsz.name) = ov;
    csAll = fr_cases27(C, kin, om);
    sel = find(strcmp({csAll.name}, gsz.name), 1);
    cs = csAll(sel);
    % signed per-constituent curves from the targets_v2 record
    T = jsondecode(fileread(fullfile(C.camp, 'targets_v2', ...
        sprintf('%s_signed_targets.json', gsz.name))));
    nG = numel(T.cells);
    % per-cell constituent table (jsondecode gives cell or struct array)
    row1 = T.cells(1);
    if iscell(row1); row1 = row1{1}; end
    mm1 = row1.muscles;
    if iscell(mm1)
        musNames = cellfun(@(x) x.muscle, mm1, 'UniformOutput', false);
    else
        musNames = arrayfun(@(x) x.muscle, mm1, 'UniformOutput', false);
    end
    nM = numel(musNames);
    tauMat = zeros(nG, nM);   % nCells x nMuscles (signed contributions)
    for cc = 1:nG
        row = T.cells(cc);
        if iscell(row); row = row{1}; end
        mm = row.muscles;
        for mi = 1:nM
            if iscell(mm)
                tauMat(cc, mi) = mm{mi}.tau_contrib_Nm;
            else
                tauMat(cc, mi) = mm(mi).tau_contrib_Nm;
            end
        end
    end

    for branch = {'rigid', 'compliant'}
        b = branch{1};
        ck = fullfile(ck1, sprintf('%s_%s.mat', gsz.name, b));
        if exist(ck, 'file') ~= 2
            fprintf(fid, '%s %s: no checkpoint - skipped\n', gsz.name, b);
            continue
        end
        S = load(ck);
        [~, rep] = fr_objective(C, cs, kin, fr_packX(cs, S.design), ...
            b, struct('quiet', true));
        tauR = rep.tauPrimary(:);
        % group-sum view (target of record)
        tauG = T.tau_vec_Nm(:);
        prjG = sign(tauG);
        prjG(abs(tauG) <= 1e-12) = 0;
        nGn = min(numel(tauR), nG);
        meetG = mean(tauR(1:nGn) .* prjG(1:nGn) >= abs(tauG(1:nGn)) - 1e-9);
        fprintf(fid, '%s %s: GROUP-sum meet %.2f (target of record)\n', ...
            gsz.name, b, meetG);
        % per named muscle
        for mi = 1:numel(musNames)
            tauM = tauMat(:, mi);
            nn = min(numel(tauR), numel(tauM));
            prjM = sign(tauM(1:nn));
            prjM(abs(tauM(1:nn)) <= 1e-12) = 0;
            meetM = mean(tauR(1:nn) .* prjM >= abs(tauM(1:nn)) - 1e-9);
            defM = max(abs(tauM(1:nn)) - abs(tauR(1:nn)));
            fprintf(fid, '  %-12s meet %.2f maxDeficit %8.3f Nm (peak |tau_m| %7.3f)\n', ...
                musNames{mi}, meetM, defM, max(abs(tauM)));
            fprintf(fcsv, '%s,%s,%s,%.4f,%.4f,%.4f,%.4f\n', ...
                gsz.name, b, musNames{mi}, meetM, defM, ...
                max(abs(tauM)), max(abs(tauR)));
        end
        nEval = nEval + 1;
    end
end
fprintf(fid, '\nevaluated %d (group, branch) pairs\n', nEval);
fprintf('per-muscle compare done: %d pairs (csv: %s)\n', nEval, csvFile);
end
