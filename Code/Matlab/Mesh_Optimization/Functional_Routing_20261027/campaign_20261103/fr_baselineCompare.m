function fr_baselineCompare(budget)
%FR_BASELINECOMPARE The 2019-vintage robotbody routings, evaluated AS-IS,
% against the signed human targets (Ben 2026-11-04: "a fair comparison of
% gait2392_robotbody's 2019 routings against human target" — the BASELINE
% the optimized rigid/compliant results must be judged against).
%
% For every group: rebuild the case with the seed geometry taken from
% gait2392_robotbody.osim's OWN path points (fr_cases27 opts.
% routeFromRobotbody), kinematics from the same robotbody model, and
% evaluate the AS-IS seed design on BOTH branches (rigid and compliant
% Xi0-Xi3). No optimization: this is the existing routing's capability.
% Scores are SIGNED, per named human muscle AND the group sum.
%
% Output: output/baseline_compare.csv (group, branch, muscle, signed
% meetFrac, maxDeficit) + output/baseline_compare.log. This is the
% baseline column of the four-way: {2019 routing baseline, optimized} x
% {rigid, compliant}, all vs the same signed human targets.
%
% Run headless: matlab -batch "addpath('<camp>'); fr_baselineCompare"
if nargin < 1
    budget = 40;   % small settle budget: not an optimization, just a
                   % taut-state settle at the as-is geometry
end

C = fr_setup();
kin = fr_kinematics(C.osimRobot);   % robotbody kinematics of record
opts = struct('routeFromRobotbody', true);
cases = fr_cases27(C, kin, containers.Map('KeyType', 'char', ...
    'ValueType', 'any'), opts);

logFile = fullfile('output', 'baseline_compare.log');
if ~exist('output', 'dir'); mkdir('output'); end
fid = fopen(logFile, 'w');
cleanupFid = onCleanup(@() fclose(fid));
fprintf(fid, ['2019 robotbody-routing baseline | %s | settle budget ' ...
    '%d evals (as-is geometry, no optimization)\n'], datetime('now'), ...
    budget);

csvFile = fullfile('output', 'baseline_compare.csv');
fcsv = fopen(csvFile, 'w');
cleanupCsv = onCleanup(@() fclose(fcsv)); %#ok<NASGU>
fprintf(fcsv, ['group,branch,muscle,signedMeetFrac,maxDeficit_Nm,' ...
    'peakAbsTauTarget_Nm,peakAbsTauBaseline_Nm\n']);

nEval = 0;
nErr = 0;
for c = 1:numel(cases)
    cs = cases(c);
    T = jsondecode(fileread(fullfile(C.camp, 'targets_v2', ...
        sprintf('%s_signed_targets.json', cs.name))));
    nG = numel(T.cells);
    row1 = T.cells(1);
    if iscell(row1); row1 = row1{1}; end
    mm1 = row1.muscles;
    if iscell(mm1)
        musNames = cellfun(@(x) x.muscle, mm1, 'UniformOutput', false);
    else
        musNames = arrayfun(@(x) x.muscle, mm1, 'UniformOutput', false);
    end
    nM = numel(musNames);
    tauMat = zeros(nG, nM);
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
    tauG = T.tau_vec_Nm(:);

    for branch = {'rigid', 'compliant'}
        b = branch{1};
        try
            so = struct('quiet', true, 'budget', budget, 'envelopes', {{}});
            [~, rep, info] = fr_solveCase(C, cs, kin, b, so);
            tauR = rep.tauPrimary(:);
            % group sum
            nn = min(numel(tauR), nG);
            prj = sign(tauG(1:nn));
            prj(abs(tauG(1:nn)) <= 1e-12) = 0;
            meetG = mean(tauR(1:nn) .* prj >= abs(tauG(1:nn)) - 1e-9);
            defG = max(abs(tauG(1:nn)) - tauR(1:nn) .* prj);
            fprintf(fcsv, '%s,%s,GROUP_SUM,%.4f,%.4f,%.4f,%.4f\n', ...
                cs.name, b, meetG, defG, max(abs(tauG(1:nn))), ...
                max(abs(tauR(1:nn))));
            fprintf(fid, '%s %s: GROUP-SUM signed meet %.2f, maxDeficit %.3f N*m\n', ...
                cs.name, b, meetG, defG);
            % per named muscle
            for mi = 1:nM
                tauM = tauMat(:, mi);
                nnM = min(numel(tauR), numel(tauM));
                prjM = sign(tauM(1:nnM));
                prjM(abs(tauM(1:nnM)) <= 1e-12) = 0;
                meetM = mean(tauR(1:nnM) .* prjM >= ...
                    abs(tauM(1:nnM)) - 1e-9);
                defM = max(abs(tauM(1:nnM)) - tauR(1:nnM) .* prjM);
                fprintf(fcsv, '%s,%s,%s,%.4f,%.4f,%.4f,%.4f\n', ...
                    cs.name, b, musNames{mi}, meetM, defM, ...
                    max(abs(tauM)), max(abs(tauR)));
            end
            nEval = nEval + 1;
        catch err
            nErr = nErr + 1;
            fprintf(fid, '%s %s: ERROR %s | %s\n', cs.name, b, ...
                err.identifier, err.message);
            fprintf(fcsv, '%s,%s,ERROR,NA,NA,NA,NA\n', cs.name, b);
        end
    end
    fprintf('baseline: %d/%d groups done (%d errors so far)\n', c, ...
        numel(cases), nErr);
end
fprintf(fid, '\nevaluated %d (group, branch) baselines, %d errors\n', ...
    nEval, nErr);
fprintf('baseline compare done: %d pairs (csv: %s)\n', nEval, csvFile);
end
