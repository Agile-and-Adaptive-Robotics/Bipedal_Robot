function fr_fourWayCompare(gridN)
%FR_FOURWAYCOMPARE Ben's 2026-11-04 four-way request: how do the BPA
% routings do with the 2019 model and with the new model, all four
% combinations (rigid/compliant x 2019/new), compared individually to
% their human specific muscle or muscle-group target.
%
% Models of record:
%   2019_masters  = Solid_Models/OpenSim/Gait2392_Robotbody/ConnorBipedal.osim
%                   (Ben's untouched master's-era robot geometry)
%   new_robotbody = Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim
%                   (the 2026-evolved robot routes copy)
% Body names are IDENTICAL between the two models (verified: same 12
% bodies), so the case machinery's frame chains resolve on both.
%
% HONEST DOMAIN NOTE (logged in the CSV too): the run_1 designs were
% OPTIMIZED with kinematics from the stock human model
% (gait2392_simbody.osim) - a documented campaign simplification. This
% comparator therefore reports out-of-sample behavior on both robot
% models; it is not a re-optimization on either.
%
% Scores are SIGNED (Ben 2026-11-04: opposite sign from the human target
% is a deficit, never a meet): per named human muscle AND per group sum.
%
% Run headless: matlab -batch "addpath('<camp>'); fr_fourWayCompare"
if nargin < 1; gridN = []; end

C = fr_setup();
szFile = fullfile(C.camp, 'sizing_winners.json');
assert(exist(szFile, 'file') == 2, 'fr_fourWayCompare:noSizing', ...
    'sizing_winners.json missing - run fr_sizingSearch first');
sz = jsondecode(fileread(szFile));
if iscell(sz.groups); groupsList = sz.groups;
else; groupsList = arrayfun(@(g) g, sz.groups, 'UniformOutput', false);
end
ck1 = fullfile(C.camp, 'run_1', 'checkpoints');

models = struct( ...
    'name', {'new_robotbody', '2019_masters'}, ...
    'osim', {fullfile(C.repo, 'Solid_Models', 'OpenSim', ...
    'Gait2392_Robotbody', 'gait2392_robotbody.osim'), ...
    fullfile(C.repo, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody', ...
    'ConnorBipedal.osim')});

logFile = fullfile('output', 'four_way_compare_v2.log');
if ~exist('output', 'dir'); mkdir('output'); end
fid = fopen(logFile, 'w');
cleanupFid = onCleanup(@() fclose(fid));
fprintf(fid, ['four-way compare | %s | gridN %d\n' ...
    'domain note: designs optimized on stock-human kinematics; both ' ...
    'robot models are out-of-sample evaluations\n'], datetime('now'), ...
    gridN);

csvFile = fullfile('output', 'four_way_compare_v2.csv');
fcsv = fopen(csvFile, 'w');
cleanupCsv = onCleanup(@() fclose(fcsv)); %#ok<NASGU>
fprintf(fcsv, ['group,branch,model,muscle,signedMeetFrac,' ...
    'maxDeficit_Nm,peakAbsTauTarget_Nm\n']);

nEval = 0;
nFail = 0;
for gi = 1:numel(groupsList)
    gsz = groupsList{gi};
    if ~isfield(gsz, 'diaStr'); continue; end
    ov = struct('diaStr', gsz.diaStr, 'nBPA', gsz.nBPA, ...
        'family', gsz.family, 'gain', gsz.gain);
    if ~isempty(gridN); ov.gridN = gridN; end
    om = containers.Map('KeyType', 'char', 'ValueType', 'any');
    om(gsz.name) = ov;

    % per-muscle signed targets from the targets_v2 record (model-
    % independent: they come from the stock human model of record)
    T = jsondecode(fileread(fullfile(C.camp, 'targets_v2', ...
        sprintf('%s_signed_targets.json', gsz.name))));
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

    for mi = 1:(nM + 1)   % last iteration = the group-sum target
        if mi <= nM
            musLabel = musNames{mi};
            tauTg = tauMat(:, mi);
        else
            musLabel = 'GROUP_SUM';
            tauTg = tauG;
        end
        for branch = {'rigid', 'compliant'}
            b = branch{1};
            ck = fullfile(ck1, sprintf('%s_%s.mat', gsz.name, b));
            if exist(ck, 'file') ~= 2
                continue
            end
            S = load(ck);
            for km = 1:numel(models)
                mname = models(km).name;
                try
                    kinM = fr_kinematics(models(km).osim);
                    casesM = fr_cases27(C, kinM, om);
                    selM = find(strcmp({casesM.name}, gsz.name), 1);
                    csM = casesM(selM);
                    [~, repM] = fr_objective(C, csM, kinM, ...
                        fr_packX(csM, S.design), b, struct('quiet', true));
                    tauR = repM.tauPrimary(:);
                    nn = min(numel(tauR), numel(tauTg));
                    prj = sign(tauTg(1:nn));
                    prj(abs(tauTg(1:nn)) <= 1e-12) = 0;
                    meetM2 = mean(tauR(1:nn) .* prj >= ...
                        abs(tauTg(1:nn)) - 1e-9);
                    defM2 = max(abs(tauTg(1:nn)) - tauR(1:nn) .* prj);
                    fprintf(fcsv, '%s,%s,%s,%s,%.4f,%.4f,%.4f\n', ...
                        gsz.name, b, mname, musLabel, meetM2, defM2, ...
                        max(abs(tauTg(1:nn))));
                    nEval = nEval + 1;
                catch err
                    nFail = nFail + 1;
                    fprintf(fid, '%s %s %s %s: ERROR %s | %s\n', ...
                        gsz.name, b, mname, musLabel, err.identifier, ...
                        err.message);
                    fprintf(fcsv, '%s,%s,%s,%s,ERROR,NA,NA\n', ...
                        gsz.name, b, mname, musLabel);
                end
            end
        end
    end
    fprintf('four-way: %s done (%d evals, %d errors so far)\n', ...
        gsz.name, nEval, nFail);
end
fprintf(fid, '\nevaluated %d model-group-branch-muscle rows, %d errors\n', ...
    nEval, nFail);
fprintf('four-way compare done: %d rows (csv: %s)\n', nEval, csvFile);
end
