% Campaign Tables 20261022
% Post-processor for the 2026-10-22 full-budget escalation-ladder campaign
% (Opt_run_BiPulley2_20261022.m, mode=escalate, OPT_BPULLEY2_RUN=1
% OPT_BPULLEY2_FULL=1). Read-only on the recorded mats/CSV; writes two NEW
% dated deliverable CSVs into Code/Matlab/Mesh_Optimization/Results:
%
%   BiPulley_escalation_choices_20261022.csv -- per-site chosen (nBPA,G)
%       with the verdict chain and the sizes solved at the chosen rung.
%   BiPulley_rom_peak_20261022.csv -- per-site RoM/peak table vs the
%       requirement envelope (graded rom_table ranges; torque target =
%       fmax x 5 cm nominal arm, the documented heuristic).
%
% Conventions carried EXACTLY from the campaign code:
%   - torque magnitude grid = vecnorm over the active crossing's 3-vector
%     (escalation_advisor_20261022.m iExtract, BiPam_pulley branch);
%   - tautMargin = min taut-feasible magnitude / tauTarget - 1 (advisor);
%   - driverMargin = the driver's st.margin recorded in the ladder CSV
%     (component-wise min over ALL finite cells, slack included -- the
%     slack-gamed number is NEVER quoted as a physical margin).
% Smoke rows (timestamps before the campaign start stamp) are excluded.

campStamp = '20261022_042213';   % campaign_start_stamp.txt
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
resultsRoot = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
ladderCsv = fullfile(resultsRoot, 'BiPulley_escalation_ladder_20261022.csv');
logFile = fullfile(scriptDir, 'bi2_fullrun_escalate_20261022.log');

T = readtable(ladderCsv, 'VariableNamingRule', 'preserve', ...
    'TextType', 'char');
ts = T.('timestamp');
% timestamps are fixed-format yyyyMMdd_HHmmss: numeric compare on the
% underscore-stripped digits (str2double on the raw string fails on '_').
campNum = str2double(strrep(campStamp, '_', ''));
keep = cellfun(@(s) str2double(strrep(s, '_', '')) >= campNum, ts);
T = T(keep, :);
fprintf('Campaign rows: %d (of %d total; smoke rows excluded).\n', ...
    height(T), numel(ts));

sites = unique(T.('name'), 'stable');
nS = numel(sites);

choicesFile = fullfile(resultsRoot, 'BiPulley_escalation_choices_20261022.csv');
romFile = fullfile(resultsRoot, 'BiPulley_rom_peak_20261022.csv');
fidC = fopen(choicesFile, 'w');
fprintf(fidC, ['site,source,rungs,configPath,chosenNBPA,chosenG,converged,' ...
    'verdictChain,driverMarginRaw,worstConstraint,tauTarget_Nm,' ...
    'minTautMag_Nm,tautMargin,peakTautMag_Nm,peakOverTarget,rest_m,' ...
    'tendon_m,kmax_m,gain,footprint_mm,budgetS,budgetP,matFile\n']);
fidR = fopen(romFile, 'w');
fprintf(fidR, ['site,joint1,joint1RangeDeg,joint2,joint2RangeDeg,' ...
    'gridCells,tautCells,slackCells,infeasibleCells,tautCoveragePct,' ...
    'joint1TautSpanDeg,joint2TautSpanDeg,failingPosesDeg,' ...
    'minTautMag_Nm,peakTautMag_Nm,tauTarget_Nm\n']);

summary = cell(nS, 10);
for s = 1:nS
    site = sites{s};
    Ts = T(strcmp(T.('name'), site), :);
    Ts = sortrows(Ts, 'rung');
    lastRow = Ts(end, :);
    % verdict chain + config path
    chain = '';
    path = '';
    for r = 1:height(Ts)
        chain = [chain, sprintf('[r%d (%d,%d) %s m %.4g t %.4g]', ...
            Ts.rung(r), Ts.nBPA(r), Ts.G(r), Ts.verdict{r}, ...
            Ts.marginRaw(r), Ts.tautMargin(r))]; %#ok<AGROW>
        if r < height(Ts)
            chain = [chain, ' -> ']; %#ok<AGROW>
        end
        if r == 1
            path = sprintf('(%d,%d)', Ts.nBPA(r), Ts.G(r));
        else
            path = sprintf('%s->(%d,%d)', path, Ts.nBPA(r), Ts.G(r));
        end
    end
    % Locate the chosen rung's mat (latest stamp if several)
    pat = sprintf('BiPulley_opensim_%s_ladder_rung%d_nB%d_G%d_*.mat', ...
        site, lastRow.rung, lastRow.nBPA, lastRow.G);
    dd = dir(fullfile(resultsRoot, pat));
    if isempty(dd)
        fprintf('MISSING MAT for %s rung %d!\n', site, lastRow.rung);
        continue;
    end
    [~, ix] = max([dd.datenum]);
    matFile = fullfile(resultsRoot, dd(ix).name);
    S = load(matFile);
    obj = S.finalObj;
    aK = obj.ActiveCrossing;
    slackMask = logical(obj.PulleySlack);
    infeasMask = logical(obj.PulleyInfeasible);
    T3 = squeeze(obj.Torque_ins(:, :, :, aK));
    mag = reshape(vecnorm(T3, 2, 2), size(slackMask));
    tautMask = ~slackMask & ~infeasMask & isfinite(mag);
    nCells = numel(slackMask);
    nSlack = nnz(slackMask);
    nInfeas = nnz(infeasMask);
    nTaut = nnz(tautMask);
    tauT = S.tauTarget;
    minTaut = min(mag(tautMask));
    peakTaut = max(mag(tautMask));
    tautMargin = minTaut / tauT - 1;
    % Requirement envelope + taut spans per joint axis
    r1 = S.spec.kin.thetaRangesDeg(1, :);
    r2 = S.spec.kin.thetaRangesDeg(2, :);
    th1 = linspace(r1(1), r1(2), size(mag, 1)) ;
    th2 = linspace(r2(1), r2(2), size(mag, 2));
    anyTaut1 = any(tautMask, 2);
    anyTaut2 = any(tautMask, 1);
    span1 = [min(th1(anyTaut1)), max(th1(anyTaut1))];
    span2 = [min(th2(anyTaut2)), max(th2(anyTaut2))];
    % failing poses (slack or infeasible), compact string
    failMask = slackMask | infeasMask;
    fp = '';
    [ii, jj] = find(failMask);
    for f = 1:numel(ii)
        kind = 'i';
        if slackMask(ii(f), jj(f)), kind = 's'; end
        fp = [fp, sprintf('%s(%.0f,%.0f)', kind, th1(ii(f)), th2(jj(f)))]; %#ok<AGROW>
        if f < numel(ii) && numel(fp) < 180, fp = [fp, ' ']; end %#ok<AGROW>
    end
    if isempty(fp), fp = 'none'; end
    if numel(fp) > 200, fp = [fp(1:197), '...']; end
    jname1 = 'knee_r'; jname2 = 'ankle_r';
    switch site
        case 'soleus_r', jname1 = 'ankle_r'; jname2 = 'fixed(0)';
        case 'tib_ant_r', jname1 = 'ankle_r'; jname2 = 'fixed(0)';
        case 'rect_fem_r', jname1 = 'hip_flex_r'; jname2 = 'fixed(0)';
        case 'bifemsh_r', jname1 = 'knee_r'; jname2 = 'fixed(0)';
    end
    fprintf(fidC, ['%s,%s,%d,"%s",%d,%d,%d,"%s",%.6g,%.6g,%.6g,%.6g,' ...
        '%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%s,%.3f,%d,%d,%s\n'], ...
        site, Ts.source{1}, height(Ts), path, lastRow.nBPA, lastRow.G, ...
        lastRow.converged, chain, lastRow.marginRaw, lastRow.worstC, ...
        tauT, minTaut, tautMargin, peakTaut, peakTaut / tauT, ...
        S.xBest(1), S.xBest(2), 0.75 * S.xBest(1), ...
        strjoin(arrayfun(@(v) sprintf('%g', v), obj.PulleyGain, ...
        'UniformOutput', false), '/'), 1000 * S.spec.footprintWidth, ...
        S.budgetS, S.budgetP, dd(ix).name);
    fprintf(fidR, ['%s,%s,[%g %g],%s,[%g %g],%d,%d,%d,%d,%.1f,' ...
        '[%.0f %.0f],[%.0f %.0f],"%s",%.6g,%.6g,%.6g\n'], ...
        site, jname1, r1, jname2, r2, nCells, nTaut, nSlack, nInfeas, ...
        100 * nTaut / nCells, span1, span2, fp, minTaut, peakTaut, tauT);
    summary(s, :) = {site, path, height(Ts), lastRow.converged, ...
        nTaut, nSlack + nInfeas, minTaut, tauT, tautMargin, dd(ix).name};
    fprintf(['%s: chosen %s after %d rung(s) (converged=%d); taut %d/%d ' ...
        'cells; min taut %.2f N*m vs target %.2f (margin %+.3f); mat %s\n'], ...
        site, path, height(Ts), lastRow.converged, nTaut, nCells, ...
        minTaut, tauT, tautMargin, dd(ix).name);
end
fclose(fidC);
fclose(fidR);
fprintf('\nWrote:\n  %s\n  %s\n', choicesFile, romFile);
