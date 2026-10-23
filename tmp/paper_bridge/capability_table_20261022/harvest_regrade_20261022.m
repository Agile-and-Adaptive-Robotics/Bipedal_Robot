% Harvest Regrade 20261022
% Assembles the dated PER-MUSCLE capability table for the 27-actuator
% survey (full run 2026-10-22 14:06-14:49, Opt_run_BiPulley_survey_
% 20261022.m) AFTER Ben's 2026-10-22 TARGET CHANGE: the human curves CSV
% was regenerated with SINGLE same-named muscle targets (no group sums).
% Per the launch ask: solves are NOT re-run; every surveyed winner is
% reloaded from its saved mat and its torque curve RE-GRADED
% deterministically against the CURRENT csv
% (human_torque_curves_20261022.csv, single-muscle, 21:59). The archived
% group-target csv (human_torque_curves_group_targets_20261022.csv,
% byte-identical to the 13:12 file the solves consumed -- same 19908
% bytes measured at 13:23) drives a verification: for every actuator
% whose target curve is UNCHANGED between the two csv versions, the
% re-grade must reproduce the checkpoint's recorded meet/active exactly;
% only changed-target actuators may move. Each row's note records which
% target version the design was OPTIMIZED against (all solves: GROUP)
% and which it is GRADED against here (SINGLE-MUSCLE).
%
% Flags: cannot_build, grade_unable, SUSPECT (meet>0 with active=0),
% capability-capped (final ladder verdict still short at the directive
% config cap {3,3}), and budget/wall cap hits (none expected -- the
% checkpoint records stopReason 'converged' for every surveyed site).
%
% PARSER NOTE: the driver-level CSVs write comma-bearing master names
% UNQUOTED ("Bicep Femoris, Long Head"), shifting those rows by one
% field; readtable chokes on the ragged rows. This harvester therefore
% uses its own quote-aware splitter (readSloppyCsv) that re-joins the
% one extra field into master_name.

scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));

ckptPath = fullfile(scriptDir, 'BiPulley_survey_full_results_20261022.csv');
ladderPath = fullfile(scriptDir, 'BiPulley_survey_full_ladder_20261022.csv');
curCsv = fullfile(scriptDir, 'human_torque_curves_20261022.csv');
grpCsv = fullfile(scriptDir, 'human_torque_curves_group_targets_20261022.csv');
for p = {ckptPath, ladderPath, curCsv, grpCsv}
    if ~isfile(p{1})
        error('harvest_regrade_20261022:missing', 'Missing input: %s', p{1})
    end
end

Ck = readSloppyCsv(ckptPath);   % struct array, one per data row
Ld = readSloppyCsv(ladderPath);
amap = jsondecode(fileread(fullfile(root, 'Code', 'Matlab', ...
    'Mesh_Optimization', 'Routing27_v2', 'actuator_map_27.json')));
acts = amap.actuators;

curvesCur = loadHumanCurves(curCsv);
curvesGrp = loadHumanCurves(grpCsv);
csvVersion = detectCsvVersion(curCsv);
fprintf('Grade CSV: %s (%s)\nOptimize-target of record: GROUP csv 13:12 (archived %s)\n', ...
    curCsv, csvVersion, grpCsv);

outFile = fullfile(scriptDir, sprintf('BiPulley_survey_perMuscle_regraded_%s.csv', ...
    char(string(datetime('now'), 'yyyyMMdd_HHmmss'))));
fid = fopen(outFile, 'w');
fprintf(fid, ['muscle,site,built,meetPct,activePct,nBPA,tackleG,' ...
    'evalsUsed,wallMin,converged,note\n']);

flagged = {};
unconverged = {};
nRows = 0;
for a = 1:numel(acts)
    A = acts(a);
    hit = find([Ck.actuator_id] == A.id, 1);
    if isempty(hit)
        error('harvest_regrade_20261022:missing', ...
            'Actuator %d (%s) missing from checkpoint', A.id, A.master_name)
    end
    R = Ck(hit);
    nRows = nRows + 1;

    if strcmp(R.status, 'cannot_build')
        fprintf(fid, '"%s",,false,0,0,0,0,0,0,false,"%s"\n', ...
            A.master_name, iCsv(['CANNOT BUILD: ' R.reason]));
        flagged{end + 1} = sprintf('%s: CANNOT BUILD', A.master_name); %#ok<AGROW>
        fprintf('%2d %-36s CANNOT BUILD\n', A.id, A.master_name)
        continue
    end
    if ~strcmp(R.status, 'surveyed')
        % grade_unable (or any other non-survey terminal status): the
        % seed BUILT but its primary dof is not modeled by the builder's
        % hinge chain -- recorded, skipped, never fatal.
        fprintf(fid, '"%s",%s,true,0,0,0,0,0,0,false,"%s"\n', ...
            A.master_name, R.seed, iCsv(['GRADE-UNABLE: ' R.reason]));
        flagged{end + 1} = sprintf('%s: GRADE-UNABLE (%s)', ...
            A.master_name, R.primary_dof); %#ok<AGROW>
        fprintf('%2d %-36s GRADE-UNABLE (%s, seed %s)\n', A.id, ...
            A.master_name, R.primary_dof, R.seed)
        continue
    end

    % ---- surveyed: deterministic re-grade of the chosen winner mat ----
    S = load(R.matFile);
    eCur = iCurvesFor(curvesCur, A.id, R.primary_dof);
    eGrp = iCurvesFor(curvesGrp, A.id, R.primary_dof);
    if isempty(eCur) || isempty(eGrp)
        error('harvest_regrade_20261022:noCurve', ...
            'No human curve for actuator %d dof %s', A.id, R.primary_dof)
    end
    targetChanged = any(abs(eCur.tauNm - eGrp.tauNm) > 1e-9);

    [meetNew, actNew] = regradeWinner(S, eCur);
    meetOld = R.meetPct;
    actOld = R.activePct;

    % Verification (the strong check): UNCHANGED target -> the re-grade
    % must reproduce the checkpoint numbers exactly (deterministic class
    % arithmetic on the saved mat, same as verify_smoke_mat_20261022).
    if ~targetChanged
        assert(abs(meetNew - meetOld) < 1e-9 && abs(actNew - actOld) < 1e-9, ...
            'harvest_regrade_20261022:mismatch', ...
            'Actuator %d: unchanged target but re-grade moved numbers.', A.id)
    end

    % Final ladder verdict (capability-capped still-short detection).
    ix = find([Ld.actuator_id] == A.id);
    [~, ord] = sort([Ld(ix).rung]);
    Ls = Ld(ix(ord));
    lastVerdict = Ls(end).verdict;

    notes = {};
    notes{end + 1} = 'optimized vs GROUP targets (full run 14:06-14:49); regraded vs SINGLE-MUSCLE csv'; %#ok<AGROW>
    if targetChanged
        notes{end + 1} = sprintf('target changed by Ben ruling: meet %.1f%%->%.1f%%, active %.1f%%->%.1f%%', ...
            meetOld, meetNew, actOld, actNew); %#ok<AGROW>
    end
    if ~strcmp(lastVerdict, 'converged')
        notes{end + 1} = sprintf('capability-capped: final rung (%d,%d) verdict %s at directive grid cap', ...
            Ls(end).nBPA, Ls(end).G, lastVerdict); %#ok<AGROW>
        flagged{end + 1} = sprintf('%s: still %s at config cap (%d,%d)', ...
            A.master_name, lastVerdict, Ls(end).nBPA, Ls(end).G); %#ok<AGROW>
    end
    if meetNew > 0 && actNew == 0
        notes{end + 1} = 'SUSPECT: meet>0 with active=0 (winner taut but zero positive-strain cells)'; %#ok<AGROW>
        flagged{end + 1} = sprintf('%s: SUSPECT meet %.1f%% with active 0%%', ...
            A.master_name, meetNew); %#ok<AGROW>
    end
    if ~R.converged || ~strcmpi(R.stopReason, 'converged')
        unconverged{end + 1} = A.master_name; %#ok<AGROW>
        notes{end + 1} = sprintf('UNCONVERGED (%s)', R.stopReason); %#ok<AGROW>
    end
    fprintf(fid, '"%s",%s,true,%.4f,%.4f,%d,%d,%d,%.2f,%d,"%s"\n', ...
        A.master_name, R.seed, meetNew, actNew, R.chosen_nB, ...
        R.chosen_G, R.evalsUsed, R.wallMin, R.converged, ...
        iCsv(strjoin(notes, '; ')));
    if targetChanged
        chgTag = ' [TARGET CHANGED]';
    else
        chgTag = '';
    end
    fprintf('%2d %-36s meet %5.1f%% (was %5.1f) active %5.1f%% (%d,%d) %s%s\n', ...
        A.id, A.master_name, meetNew, meetOld, actNew, R.chosen_nB, ...
        R.chosen_G, lastVerdict, chgTag)
end
fclose(fid);

assert(nRows == 27, 'harvest_regrade_20261022:rowCount', ...
    'Wrote %d rows, expected 27', nRows)
fprintf('\n==== PER-MUSCLE CAPABILITY TABLE (regraded, 27 rows) ====\n');
type(outFile);
fprintf('\nRows: %d. Flags (%d):\n', nRows, numel(flagged));
for f = 1:numel(flagged)
    fprintf('  - %s\n', flagged{f});
end
fprintf('Budget/wall-cap unconverged (%d): %s\n', numel(unconverged), ...
    strjoin(unconverged, ', '));
fprintf('WROTE %s\n', outFile);

%% ---------------------------------------------------------------------
function rows = readSloppyCsv(path)
% Quote-aware comma splitter for the driver-level CSVs. The driver wrote
% master names like "Bicep Femoris, Long Head" UNQUOTED, so those rows
% carry one extra field; the extra field is re-joined into the column
% named master_name (the only comma-bearing unquoted column).
fid = fopen(path, 'r');
lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
lines = lines{1};
hdr = splitCsvLine(lines{1});
nCol = numel(hdr);
cols = cellfun(@(c) strtrim(stripq(c)), hdr, 'UniformOutput', false);
rows = struct([]);
n = 0;
for i = 2:numel(lines)
    if isempty(strtrim(lines{i}))
        continue
    end
    f = splitCsvLine(lines{i});
    if numel(f) == nCol + 1 && any(strcmpi(cols, 'master_name'))
        mi = find(strcmpi(cols, 'master_name'), 1);
        f = [f(1:mi - 1); {sprintf('%s,%s', f{mi}, f{mi + 1})}; f(mi + 2:end)];
    end
    if numel(f) ~= nCol
        error('harvest_regrade_20261022:parse', ...
            'Row %d of %s has %d fields (header %d) after the merge rule.', ...
            i, path, numel(f), nCol)
    end
    n = n + 1;
    if n == 1
        rows = struct();   % field-assignment into struct([]) errors
    end
    for c = 1:nCol
        v = strtrim(stripq(f{c}));
        if ismember(lower(cols{c}), {'actuator_id', 'rung', 'nbpa', 'g', ...
                'chosen_nb', 'chosen_g', 'evalsused', 'slackcells', ...
                'infeasiblecells', 'tautcells'})
            rows(n).(fixName(cols{c})) = str2double(v); %#ok<AGROW>
        elseif ismember(lower(cols{c}), {'meetpct', 'activepct', 'score', ...
                'wallmin', 'rest_m', 'tendon_m', 'mintauttau_nm', ...
                'targettau_nm'})
            rows(n).(fixName(cols{c})) = str2double(v); %#ok<AGROW>
        elseif ismember(lower(cols{c}), {'converged'})
            rows(n).(fixName(cols{c})) = str2double(v) ~= 0; %#ok<AGROW>
        else
            rows(n).(fixName(cols{c})) = v; %#ok<AGROW>
        end
    end
end
end

function f = splitCsvLine(line)
% Split on commas outside double quotes; keep quote contents intact.
f = {};
cur = '';
inQ = false;
for i = 1:numel(line)
    ch = line(i);
    if ch == '"'
        inQ = ~inQ;
        cur = [cur, ch]; %#ok<AGROW>
    elseif ch == ',' && ~inQ
        f{end + 1, 1} = cur; %#ok<AGROW>
        cur = '';
    else
        cur = [cur, ch]; %#ok<AGROW>
    end
end
f{end + 1, 1} = cur; %#ok<AGROW>
end

function s = stripq(s)
if numel(s) >= 2 && s(1) == '"' && s(end) == '"'
    s = s(2:end - 1);
end
end

function nm = fixName(c)
nm = regexprep(regexprep(c, '[^A-Za-z0-9_]', '_'), '_+', '_');
if isempty(nm) || ~isletter(nm(1))
    nm = ['x', nm];
end
end

function [meetPct, activePct] = regradeWinner(S, e)
% Deterministic re-grade of a saved winner: rebuild the class at the
% mat's own xBest/T/spec and grade against curve e on the mat's grid
% (the exact runRung/grade arithmetic of the driver, lines 502-532).
obj = BiPam_pulley2_20261022(S.spec.name, S.points, S.spec.cross, 20, ...
    S.T, S.xBest(1), 0.75 * S.xBest(1), S.xBest(2), 0.025, 620, ...
    0.005, 5e5, 1e4, 6, S.spec.pulley);
slot = S.surveyMeta.slot;
axisU = S.spec.kin.axes(slot, :) / norm(S.spec.kin.axes(slot, :));
[N1, N2] = obj.gridSize();
T3 = obj.Torque_ins(:, :, :, slot);
tauB = reshape(sum(T3 .* reshape(axisU, 1, 3, 1), 2), N1, N2);
strain = reshape(obj.strain_p, N1, N2);
slack = logical(obj.PulleySlack);
infeas = logical(obj.PulleyInfeasible);
taut = ~slack & ~infeas & isfinite(tauB);
% Target-grid orientation EXACTLY as the full-run driver (lines 388-397):
% slot-1 swept -> primary pose along dim 1 (N2 = duplicated pinned cols);
% slot-2 swept -> N1 = N2 = 25 with the pinned joint-1 transform
% duplicated across dim 3, primary pose along dim 2 (columns).
if slot == 1
    tauH = repmat(reshape(e.tauNm, [], 1), 1, N2);
else
    tauH = repmat(reshape(e.tauNm, 1, []), N1, 1);
end
meet = taut & isfinite(tauH) & ...
    (sign(tauB) == sign(tauH) | abs(tauH) <= 1e-9) & ...
    abs(tauB) >= abs(tauH) - 1e-12;
active = taut & isfinite(strain) & strain > 0 & abs(tauB) > 0;
meetPct = 100 * sum(meet(:)) / numel(meet);
activePct = 100 * sum(active(:)) / numel(active);
end

function v = detectCsvVersion(csvPath)
% Version tag from the csv's leading '#' comments (the regenerated file
% announces SINGLE-MUSCLE; the archived group file says human_group).
fid = fopen(csvPath, 'r');
txt = fread(fid, '*char')';
fclose(fid);
hdr = txt(1:min(numel(txt), 2000));
if contains(hdr, 'SINGLE-MUSCLE', 'IgnoreCase', true)
    v = 'single-muscle';
elseif contains(hdr, 'human_group')
    v = 'group';
else
    v = 'unknown';
end
end

function curves = loadHumanCurves(csvPath)
% Same wide loader as the driver (header inspected 2026-10-22): one row
% per actuator -- actuator_id, master_name, primary_dof,
% theta_deg_01..25, tau_h_Nm_01..25, artifact_flag_caveat; '#' comments.
raw = readtable(csvPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
cols = raw.Properties.VariableNames;
angCols = iOrderedCols(cols, 'theta_deg_');
tauCols = iOrderedCols(cols, 'tau_h_Nm_');
curves = containers.Map('KeyType', 'double', 'ValueType', 'any');
for r = 1:height(raw)
    idn = double(raw.actuator_id(r));
    angles = zeros(1, numel(angCols));
    taus = zeros(1, numel(tauCols));
    for c = 1:numel(angCols)
        angles(c) = raw.(angCols{c})(r);
        taus(c) = raw.(tauCols{c})(r);
    end
    e = struct('dof', strtrim(raw.primary_dof{r}), ...
        'anglesDeg', angles, 'tauNm', taus);
    if isKey(curves, idn)
        curves(idn) = [curves(idn), e];
    else
        curves(idn) = e;
    end
end
end

function picked = iOrderedCols(cols, prefix)
idx = [];
for c = 1:numel(cols)
    if strncmpi(cols{c}, prefix, numel(prefix))
        n = str2double(strtrim(cols{c}(numel(prefix) + 1:end)));
        if ~isnan(n)
            idx(end + 1, :) = [n, c];   %#ok<AGROW>
        end
    end
end
idx = sortrows(idx, 1);
picked = cols(idx(:, 2));
end

function e = iCurvesFor(curves, actId, dof)
e = [];
if isKey(curves, double(actId))
    es = curves(double(actId));
else
    return
end
for k = 1:numel(es)
    if strcmp(es(k).dof, dof)
        e = es(k);
        return
    end
end
end

function s = iCsv(txtIn)
s = strrep(txtIn, '"', '''');
end
