% study2_capability_v3_20261027.m -- CAPABILITY COMPARISON v3 (2026-10-27)
% Three-metric masters-vs-updated assembly across ALL 27 actuators, per
% the 2026-10-27 studies directive (the grading script study2's header
% said was never written; canonical function-file layout like study2).
%
% SIDES:
%   UPDATED  = the survey2 LADDER WINNER of record per site
%              (BiPulley_survey2_perMuscle_20261023.csv matFile column;
%              masters routing, old-code diameter, config escalated
%              (nB,G) with sizes re-solved per rung).
%   ORIGINAL = the masters design intent at its base configuration:
%              - the EIGHT 40 mm sites: STUDY 2 outputs (survey3_orig
%                mats: masters route verbatim, the masters 40 mm BPA
%                substituted by FOUR PARALLEL 20 mm BPAs, survey2
%                rung-1 protocol) -- BiPulley_survey3_orig_4x20_20261027.csv;
%              - the OTHER 19 sites: the SAME SURVEY2 RUN AT RUNG 1
%                (nB1, G1 mats: masters routing + old-code diameter,
%                sizes re-solved at rung 1). NOTE, stated plainly: these
%                are the survey2 ladder's own rung-1 mats -- the same
%                masters-routing machinery as the updated side; wherever
%                the ladder winner IS rung 1 the very same mat would
%                serve both sides (as it happens no 19-site winner is
%                rung 1 this run; the flag same_mat records it per row).
%
% METRICS (per actuator, on the slot-layout 5x5 graded grid, DOMAIN =
% in_graded_domain mask from human_surfaces_5x5_grading_20261023.csv as
% carried in each mat's dom; domain-excluded cells are OUTSIDE every
% denominator per the grading contract, never counted as failures):
%   M1 boolean area share = 100 * sum(meet & dom) / sum(dom), meet =
%     the drivers' exact boolean (taut & sign-agree & |tauR| >= |tauH|).
%   M2 average of the two axis means = 100 * (mean_i rowFrac_i +
%     mean_j colFrac_j) / 2, rowFrac_i = meet share within row i's
%     domain cells (rows/cols with no domain cell omitted from that
%     axis mean). On FULL rectangular domains M2 == M1 identically
%     (Fubini); they differ only where the domain is clipped (the
%     hip-envelope sites) -- asserted numerically below.
%   M3 mean margin = mean over DOMAIN cells of min(tauR/tauH, 1),
%     SIGN-AWARE (sign disagreement gives a negative ratio; slack cells
%     give 0; cap 1 so excess never inflates), with dead cells COUNTED:
%     cells whose robot torque is NaN (infeasible) enter as 0 rather
%     than being dropped from the denominator.
%
% Moment-arm rules are ALREADY applied inside the mats' graded surfaces
% (both sides graded against the same human_surfaces CSV); this script
% only reuses them.
%
% Outputs (dated 20261027):
%   tmp/paper_bridge/capability_table_20261022/capability_v3_three_metric_20261027.csv
%   Results/BiPulley_surveys/Capability_v3_three_metric_20261027.mat

function study2_capability_v3_20261027()
R = 'D:\GitHub\Bipedal_Robot';
capDir = fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022');
resDir = fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_surveys');
addpath(genpath(fullfile(R, 'Code', 'Matlab')));

%% Roster: survey2 winners of record -------------------------------------
Tpm = readtable(fullfile(capDir, 'BiPulley_survey2_perMuscle_20261023.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
assert(height(Tpm) == 27, 'perMuscle roster %d rows', height(Tpm));
assert(all(strcmp(Tpm.built, 'true')), 'a site is built:false');

%% The eight substitution sites (Study 2) --------------------------------
T3 = readtable(fullfile(capDir, 'BiPulley_survey3_orig_4x20_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
subSite = T3.site;
subMat = T3.matFile;
assert(numel(subSite) == 8, 'survey3 roster %d rows', numel(subSite));

%% Sweep all 27 -----------------------------------------------------------
hdr = {'actuator_id', 'site', 'orig_source', 'orig_nBPA', 'orig_G', ...
    'orig_dia_mm', 'upd_nBPA', 'upd_G', 'upd_dia_mm', ...
    'orig_M1_areashare', 'orig_M2_axismean', 'orig_M3_margin', ...
    'upd_M1_areashare', 'upd_M2_axismean', 'upd_M3_margin', ...
    'dM1', 'dM2', 'dM3', 'orig_meetCells', 'orig_domCells', ...
    'upd_meetCells', 'upd_domCells', 'fubini_absdiff_orig', ...
    'fubini_absdiff_upd', 'same_mat'};
rows = cell(27, numel(hdr));
fubMax = struct('orig', 0, 'upd', 0);
progFid = fopen(fullfile(capDir, 'capability_v3_progress_20261027.log'), 'w');
tAll = tic;
for k = 1:27
    site = strtrim(Tpm.site{k});
    idn = double(Tpm.actuator_id(k));
    % Winner-mat path repair: the perMuscle CSV was written BEFORE the
    % 10-23 BiPulley_surveys consolidation, so its matFile column points
    % at the retired folder ...\Results\BiPulley_survey\<file>. Resolve
    % by BASENAME in the current consolidated folder (loud if absent).
    updWant = strtrim(Tpm.matFile{k});
    [~, bn, bx] = fileparts(updWant);
    updMat = fullfile(resDir, [bn, bx]);
    if ~strcmp(updWant, updMat) && isfile(updWant)
        updMat = updWant;               % pre-consolidation path still live
    end
    assert(isfile(updMat), '%s: winner mat %s (CSV %s) not found', ...
        site, updMat, updWant);

    ixSub = find(strcmp(subSite, site), 1);
    if ~isempty(ixSub)
        origMat = strtrim(subMat{ixSub});
        origTag = 'survey3_orig_4x20';
    else
        d = dir(fullfile(resDir, ['BiPulley_survey2_', site, ...
            '_nB1_G1_*.mat']));
        assert(~isempty(d), '%s: no survey2 rung-1 mat', site);
        [~, im] = max([d.datenum]);
        origMat = fullfile(d(im).folder, d(im).name);
        origTag = 'survey2_rung1';
    end

    Mo = loadMatMetrics(origMat);
    Mu = loadMatMetrics(updMat);
    fprintf(progFid, ['[%d] %-12s orig %-16s meet %2d/%2d (rec %g) | ' ...
        'upd meet %2d/%2d (rec %g)\n'], idn, site, origTag, ...
        Mo.meetCells, Mo.domCells, Mo.recMeet, Mu.meetCells, ...
        Mu.domCells, Mu.recMeet);
    fclose(progFid);
    progFid = fopen(fullfile(capDir, 'capability_v3_progress_20261027.log'), 'a');

    % Fubini identity check on full rectangular domains
    fo = NaN; fu = NaN;
    if Mo.domCells == 25
        fo = abs(Mo.M1 - Mo.M2);
        assert(fo < 1e-9, '%s orig: full domain but |M1-M2|=%.3g', site, fo);
        fubMax.orig = max(fubMax.orig, fo);
    end
    if Mu.domCells == 25
        fu = abs(Mu.M1 - Mu.M2);
        assert(fu < 1e-9, '%s upd: full domain but |M1-M2|=%.3g', site, fu);
        fubMax.upd = max(fubMax.upd, fu);
    end

    rows(k, :) = {idn, site, origTag, Mo.nB, Mo.G, Mo.dia, Mu.nB, Mu.G, ...
        Mu.dia, Mo.M1, Mo.M2, Mo.M3, Mu.M1, Mu.M2, Mu.M3, ...
        Mu.M1 - Mo.M1, Mu.M2 - Mo.M2, Mu.M3 - Mo.M3, ...
        Mo.meetCells, Mo.domCells, Mu.meetCells, Mu.domCells, fo, fu, ...
        strcmp(origMat, updMat)};
    fprintf(['[%2d] %-12s ORIG(%s %dx%d d%dm) M1 %6.2f M2 %6.2f M3 %+.4f' ...
        ' | UPD(%dx%d d%dm) M1 %6.2f M2 %6.2f M3 %+.4f\n'], idn, site, ...
        origTag, Mo.nB, Mo.G, Mo.dia, Mo.M1, Mo.M2, Mo.M3, ...
        Mu.nB, Mu.G, Mu.dia, Mu.M1, Mu.M2, Mu.M3);
end
fprintf('Assembly wall time: %.1f s\n', toc(tAll));
fclose(progFid);

%% Row averages + summary -------------------------------------------------
M1o = mean([rows{:, 10}]); M2o = mean([rows{:, 11}]); M3o = mean([rows{:, 12}]);
M1u = mean([rows{:, 13}]); M2u = mean([rows{:, 14}]); M3u = mean([rows{:, 15}]);
nUp1 = nnz(cell2mat(rows(:, 16)) > 1e-9);
nDn1 = nnz(cell2mat(rows(:, 16)) < -1e-9);
nUp3 = nnz(cell2mat(rows(:, 18)) > 1e-9);
nDn3 = nnz(cell2mat(rows(:, 18)) < -1e-9);
nonRect = find(cell2mat(rows(:, 20)) < 25 | cell2mat(rows(:, 22)) < 25);
m12max = 0;
for k = nonRect(:).'
    m12max = max([m12max, abs([rows{k, 10}] - [rows{k, 11}]), ...
        abs([rows{k, 13}] - [rows{k, 14}])]);
end
fprintf('\n================ CAPABILITY v3 THREE-METRIC SUMMARY ================\n');
fprintf('row averages over 27 actuators:\n');
fprintf('  ORIGINAL : M1 %7.3f%%  M2 %7.3f%%  M3 %+.4f\n', M1o, M2o, M3o);
fprintf('  UPDATED  : M1 %7.3f%%  M2 %7.3f%%  M3 %+.4f\n', M1u, M2u, M3u);
fprintf('delta    : dM1 %+.3f  dM2 %+.3f  dM3 %+.4f\n', ...
    M1u - M1o, M2u - M2o, M3u - M3o);
fprintf(['per-metric updated-vs-original wins: M1 up %d / down %d ; ' ...
    'M3 up %d / down %d (of 27)\n'], nUp1, nDn1, nUp3, nDn3);
fprintf(['Fubini max |M1-M2| on full-domain cells: orig %.2g, upd %.2g ' ...
    '(identity holds)\n'], fubMax.orig, fubMax.upd);
fprintf(['non-rectangular domains (M1 != M2 possible): %d site(s); max ' ...
    '|M1-M2| there %.3f pct-pts\n'], numel(nonRect), m12max);

%% CSV ---------------------------------------------------------------------
outCsv = fullfile(capDir, 'capability_v3_three_metric_20261027.csv');
fid = fopen(outCsv, 'w');
fprintf(fid, ['# capability_v3_three_metric_20261027 - masters-vs-updated ' ...
    'three-metric comparison, all 27 actuators (2026-10-27)\n']);
fprintf(fid, ['# UPDATED = survey2 ladder winners (masters routing, ' ...
    'old-code dia, escalated config); ORIGINAL = 8 x 40mm sites from ' ...
    'study2 survey3_orig mats (4x 20mm substitution, masters route)\n']);
fprintf(fid, ['# + 19 other sites from their survey2 RUNG-1 (nB1,G1) mats ' ...
    '(masters routing + old-code diameters; the same mats that serve the ' ...
    'updated side at rung-1 -- same_mat flags any row where the winner ' ...
    'mat IS the rung-1 mat)\n']);
fprintf(fid, ['# M1 boolean area share pct; M2 average of the two axis ' ...
    'means pct; M3 mean margin = mean over domain cells of ' ...
    'min(tauR/tauH,1), sign-aware, dead (NaN-torque) cells counted as 0\n']);
fprintf(fid, ['# domain mask per the grading surfaces CSV ' ...
    '(human_surfaces_5x5_grading_20261023.csv) as carried in each mat; ' ...
    'domain-excluded cells outside every denominator\n']);
fprintf(fid, '%s\n', strjoin(hdr, ','));
for k = 1:27
    fprintf(fid, ['%d,%s,%s,%d,%d,%d,%d,%d,%d,%.4f,%.4f,%.6f,%.4f,%.4f,' ...
        '%.6f,%+.4f,%+.4f,%+.6f,%d,%d,%d,%d,%s,%s,%d\n'], ...
        rows{k, 1}, rows{k, 2}, rows{k, 3}, rows{k, 4}, rows{k, 5}, ...
        rows{k, 6}, rows{k, 7}, rows{k, 8}, rows{k, 9}, rows{k, 10}, ...
        rows{k, 11}, rows{k, 12}, rows{k, 13}, rows{k, 14}, rows{k, 15}, ...
        rows{k, 16}, rows{k, 17}, rows{k, 18}, rows{k, 19}, rows{k, 20}, ...
        rows{k, 21}, rows{k, 22}, fmtOrNa(rows{k, 23}), ...
        fmtOrNa(rows{k, 24}), rows{k, 25});
end
fprintf(fid, '# row averages,ORIGINAL,%.4f,%.4f,%.6f\n', M1o, M2o, M3o);
fprintf(fid, '# row averages,UPDATED,%.4f,%.4f,%.6f\n', M1u, M2u, M3u);
fprintf(fid, ['# row averages,DELTA,%+.4f,%+.4f,%+.6f ; updated wins M1 ' ...
    '%d/%d, M3 %d/%d ; Fubini max|M1-M2| full domains orig %.2g upd ' ...
    '%.2g\n'], M1u - M1o, M2u - M2o, M3u - M3o, nUp1, 27, nUp3, 27, ...
    fubMax.orig, fubMax.upd);
fclose(fid);
fprintf('CSV written: %s\n', outCsv);

save(fullfile(resDir, 'Capability_v3_three_metric_20261027.mat'), ...
    'rows', 'hdr', 'M1o', 'M2o', 'M3o', 'M1u', 'M2u', 'M3u', 'fubMax');
fprintf('CAPABILITY v3 DONE\n');
end

%% --------------- one mat -> the three metrics --------------------------
function M = loadMatMetrics(matFile)
S = load(matFile, 'surveyMeta', 'tauHuman', 'dom', 'g');
assert(isfield(S, 'g') && isfield(S.g, 'tauB'), '%s: no g.tauB', matFile);
tauH = S.tauHuman;
dom = logical(S.dom);
tauR = S.g.tauB;
taut = ~logical(S.g.slack) & ~logical(S.g.infeas) & isfinite(tauR);
meet = taut & isfinite(tauH) & ...
    (sign(tauR) == sign(tauH) | abs(tauH) <= 1e-9) & ...
    abs(tauR) >= abs(tauH) - 1e-12;
nDom = sum(dom(:));
meetCells = sum(meet(:) & dom(:));   % (:) repair: sum of a 5x5 logical is a
                                     % column vector, not a total (the
                                     % assert-cond mystery, 2026-10-27)

% Cross-check against the mat's own recorded count (both mat families
% store surveyMeta.meetCells at xBest). Precise, loud diagnostics: a
% non-scalar/non-finite record or a genuine mismatch each aborts with
% the mat named (nothing silent).
M.recMeet = NaN;
if isfield(S.surveyMeta, 'meetCells')
    rec = S.surveyMeta.meetCells;
    assert(isscalar(rec) && isfinite(rec) && isnumeric(rec), ...
        '%s: meetCells record not scalar-finite numeric (%s %s)', ...
        matFile, class(rec), mat2str(size(rec)));
    M.recMeet = double(rec);
    assert(meetCells == M.recMeet, ...
        '%s: recomputed meet %d vs recorded %g', matFile, meetCells, ...
        M.recMeet);
end

M1 = 100 * meetCells / nDom;

% M2: per-axis meet fractions within each row/col's domain cells.
[nR, nC] = size(dom);
fr = []; fc = [];
for i = 1:nR
    c = sum(dom(i, :));
    if c > 0
        fr(end + 1) = sum(meet(i, :) & dom(i, :)) / c; %#ok<AGROW>
    end
end
for j = 1:nC
    c = sum(dom(:, j));
    if c > 0
        fc(end + 1) = sum(meet(:, j) & dom(:, j)) / c; %#ok<AGROW>
    end
end
M2 = 100 * (mean(fr) + mean(fc)) / 2;

% M3: sign-aware capped ratio; dead (NaN) cells counted as 0.
v = min(tauR ./ tauH, 1);
v(~isfinite(v)) = 0;
M3 = mean(v(dom));

if isfield(S.surveyMeta, 'bpa_size_mm')
    dia = double(S.surveyMeta.bpa_size_mm);
elseif isfield(S.surveyMeta, 'dia_mm')
    dia = double(S.surveyMeta.dia_mm);
else
    dia = NaN;
end
M = struct('M1', M1, 'M2', M2, 'M3', M3, 'meetCells', meetCells, ...
    'domCells', nDom, 'nB', double(S.surveyMeta.nBPA), ...
    'G', double(S.surveyMeta.G), 'dia', dia, 'recMeet', M.recMeet);
end

function s = fmtOrNa(x)
if isfinite(x)
    s = sprintf('%.2g', x);
else
    s = 'NA';
end
end
