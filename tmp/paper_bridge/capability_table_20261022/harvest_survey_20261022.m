% Harvest Survey 20261022
% Assembles the dated PER-MUSCLE capability table (one row per map
% actuator, built:false rows for muscles that could not build) from the
% full-run checkpoint CSV written by Opt_run_BiPulley_survey_20261022.m
% (BiPulley_survey_full_results_20261022.csv), joined against
% actuator_map_27.json for the map order. Flags: cannot-build,
% grade-unable, cap-hit-unconverged, and suspicious rows
% (meetPct > 0 while activePct = 0; all-dead 0%/0% surveys). Verifies
% 27 rows, one per actuator, none missing. Writes
% BiPulley_survey_perMuscle_<stamp>.csv next to this script and prints
% the full table + flag summary.
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    root = parent;
end

ckptPath = fullfile(scriptDir, 'BiPulley_survey_full_results_20261022.csv');
if ~isfile(ckptPath)
    error('harvest_survey_20261022:noCkpt', ...
        'Checkpoint CSV not found: %s', ckptPath)
end
T = readtable(ckptPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char');
amap = jsondecode(fileread(fullfile(root, 'Code', 'Matlab', ...
    'Mesh_Optimization', 'Routing27_v2', 'actuator_map_27.json')));
acts = amap.actuators;
nA = numel(acts);

outFile = fullfile(scriptDir, sprintf('BiPulley_survey_perMuscle_%s.csv', ...
    char(string(datetime('now'), 'yyyyMMdd_HHmmss'))));
fid = fopen(outFile, 'w');
fprintf(fid, ['muscle,site,built,meetPct,activePct,nBPA,tackleG,' ...
    'evalsUsed,wallMin,converged,note\n']);

flagged = {};
unconverged = {};
nMissing = 0;
for a = 1:nA
    A = acts(a);
    hit = find(double(T.actuator_id) == A.id, 1);
    muscle = A.master_name;
    if isempty(hit)
        nMissing = nMissing + 1;
        fprintf(fid, '%s,,false,,,,,,,,MISSING FROM CHECKPOINT\n', muscle);
        flagged{end + 1} = sprintf('%s: MISSING from checkpoint', muscle); %#ok<AGROW>
        unconverged{end + 1} = muscle; %#ok<AGROW>
        continue
    end
    R = T(hit, :);
    status = strtrim(R.status);
    site = strtrim(R.seed);
    reason = strtrim(R.reason);
    stopReason = strtrim(R.stopReason);
    if strcmp(status, 'surveyed')
        built = true;
        meetPct = R.meetPct;
        activePct = R.activePct;
        nB = R.chosen_nB;
        G = R.chosen_G;
        evals = R.evalsUsed;
        wall = R.wallMin;
        conv = logical(R.converged);
        notes = {};
        if ~conv
            notes{end + 1} = sprintf('UNCONVERGED (%s)', stopReason); %#ok<AGROW>
            unconverged{end + 1} = muscle; %#ok<AGROW>
        end
        if meetPct > 0 && activePct == 0
            notes{end + 1} = 'SUSPECT: meet>0 with active=0'; %#ok<AGROW>
            flagged{end + 1} = sprintf('%s: SUSPECT meet %.1f%% with active 0%%', ...
                muscle, meetPct); %#ok<AGROW>
        end
        if activePct == 0
            notes{end + 1} = sprintf('all-dead survey (%s)', stopReason); %#ok<AGROW>
            flagged{end + 1} = sprintf( ...
                '%s: 0%% active (no taut cell in bounds on the graded grid)', ...
                muscle); %#ok<AGROW>
        end
        note = strjoin(notes, '; ');
        fprintf(fid, ['%s,%s,true,%.4f,%.4f,%d,%d,%d,%.2f,%d,"%s"\n'], ...
            muscle, site, meetPct, activePct, nB, G, evals, wall, conv, note);
    elseif strcmp(status, 'cannot_build')
        built = false;
        fprintf(fid, '%s,,false,,,,,,,,CANNOT BUILD: "%s"\n', muscle, reason);
        flagged{end + 1} = sprintf('%s: CANNOT BUILD (%s)', muscle, reason); %#ok<AGROW>
    else
        % grade_unable | no_human_curve | anything else: built but not
        % surveyable -- recorded with the reason.
        fprintf(fid, ['%s,%s,true,,,,,,,,%s: "%s"\n'], muscle, site, ...
            upper(status), reason);
        flagged{end + 1} = sprintf('%s: %s (%s)', muscle, status, reason); %#ok<AGROW>
    end
end
fclose(fid);

fprintf('\n==== PER-MUSCLE CAPABILITY TABLE (27 actuators) ====\n');
type(outFile);
fprintf('\nRow count check: %d map actuators, %d checkpoint rows, %d missing.\n', ...
    nA, height(T), nMissing);
fprintf('Flagged (%d):\n', numel(flagged));
for f = 1:numel(flagged)
    fprintf('  - %s\n', flagged{f});
end
fprintf('Unconverged (%d): %s\n', numel(unconverged), ...
    strjoin(unconverged, ', '));
fprintf('WROTE %s\n', outFile);
if nMissing > 0
    error('harvest_survey_20261022:incomplete', ...
        '%d actuator(s) missing from the checkpoint -- rerun the driver.', nMissing)
end
