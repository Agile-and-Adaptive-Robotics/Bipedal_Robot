function [advice, report] = escalation_advisor_20261022(objOrMat, targetTau, opts)
% Escalation Advisor 20261022
% Fuzzy configuration-escalation advisor per the pulley sweep directive
% section 5 (Ben Bolen, 2026-10-21: "Lacking in torque? -> more BPAs in
% parallel. Lacking in range of motion? -> more tackle line parts."), per
% site and per pose region -- OPTIMIZER_IMPROVEMENT_PLAN.md Item 2 with
% the 2026-10-22 review blocker fix applied.
%
%   BLOCKER FIX (review 2026-10-22): the torque-short test is defined
%   over TAUT, FEASIBLE cells ONLY. Torque_ins is identically 0 at slack
%   cells by construction (F_ins zeroed at slack: MonoPam_pulley.m:595
%   "F_ins_i(slack_i, :) = 0"; BiPam_pulley.m:671 maskZero3 on
%   PulleySlack), so a test that includes slack cells reads ANY
%   slack-bearing run as torque-short regardless of its taut-cell torque.
%   The plan's own replay fixture proved it:
%   doctrine_rerun_20261021.log:1160 records margin = -1.0000 for
%   med_gas_r at (2,[2 1]) -- margin -1 means min(validTau) = 0, i.e. the
%   binding cell is a SLACK cell -- with 'slack cells = 10' and worst
%   constraint +0.430868 in the same log block. Under this advisor the
%   raw (slack-gamed) margin -1 is REPORTED but never consumed; the
%   verdict's torque side fires only on genuine taut-cell shortfalls.
%   That fixture's honest verdict under the fixed rule is BOTH (10 slack
%   cells AND 5 genuinely short taut cells at 7.6-28.8 N*m vs the
%   77.9 N*m target), i.e. it steps BOTH sides -- critically NOT the
%   plan's slack-driven verdict 'torque-short -> advise nBPA 3 alone',
%   and a taut-clean design would read rom_slack_short alone (unit test
%   (a)). SLACK and INFEASIBLE cells always map to the RoM/slack side.
%
% INPUTS
%   objOrMat   one of:
%     - BiPam_pulley object: torque = |Torque_ins(:,:,ActiveCrossing)|,
%       masks PulleySlack / PulleyInfeasible (N1 x N2);
%     - MonoPam_pulley object: torque = |Torque_ins(:,3)| (the flexor
%       pipeline's TorqueInsZ convention, Opt_run_pulley.m:1200-1206),
%       masks PulleySlack / PulleyInfeasible (N x 1);
%     - struct with fields TorqueIns (grid), PulleySlack, PulleyInfeasible
%       (logical grids, same size) and optional ActiveCrossing -- the
%       synthetic-test entry (no class construction needed);
%     - char mat path of an Opt_run_BiPulley result mat: loads finalObj
%       (+ tauTarget when targetTau is omitted).
%   targetTau  scalar torque target (N*m). [] -> taken from the mat's
%      tauTarget field (mat input only).
%   opts       struct with optional fields:
%      .verbose     (default true) print the diagnosis block
%      .capNBPA     (default 3) directive grid cap nBPA <= 3
%      .capG        (default 3) directive grid cap G <= 3
%      .curNBPA     (default 1) current ladder position (stepping advice)
%      .curG        (default 1)
%      .cellWeights W grid (default ones): cells with W <= 0 are excluded
%      from BOTH tests (the Item 6 soft-corner mask hook: a ladder step
%      cannot be justified by an artifact-corner cell).
%      .csvFile     append one report row to this CSV (created w/ header)
%
% OUTPUTS
%   advice  struct:
%     .verdict    'converged' | 'torque_short' | 'rom_slack_short' | 'both'
%     .stepNBPA   0/1  (+1 = add one parallel BPA; 0 at grid cap)
%     .stepG      0/1  (+1 = add one tackle line part; 0 at grid cap)
%     .nextNBPA   nBPA after the step (NaN if stepNBPA == 0)
%     .nextG      G after the step (NaN if stepG == 0)
%     .converged  logical
%     .reason     char chain justifying the choice (the logged choice)
%   report  struct (anti-gaming: raw AND taut margins plus the
%     per-region breakdown, so no step can hide behind a slack-artifact
%     cell): nCells, nSlack, nInfeasible, nTautFeasible, minTautTau,
%     tautMargin, rawMargin (the slack-gamed min over ALL finite cells,
%     reported to expose the artifact), minTautCell ([i j] or [i]),
%     torqueTestEvaluable, worstSlackSide ('slack'/'infeasible'/''/both),
%     targetTau, source.
%
% Deterministic: zero optimizer calls, zero RNG. Read-only on any input.

% --- Input resolution ---------------------------------------------------
if nargin < 2, targetTau = []; end
if nargin < 3, opts = struct(); end
verbose = ~isfield(opts, 'verbose') || opts.verbose;
capNBPA = iGet(opts, 'capNBPA', 3);
capG    = iGet(opts, 'capG', 3);
curNBPA = iGet(opts, 'curNBPA', 1);
curG    = iGet(opts, 'curG', 1);

srcLabel = '';
if ischar(objOrMat) || isstring(objOrMat)
    srcLabel = char(objOrMat);
    S = load(srcLabel);
    if isempty(targetTau) && isfield(S, 'tauTarget')
        targetTau = S.tauTarget;
    end
    objOrMat = S.finalObj;
end

[TorqueIns, slackMask, infeasMask] = iExtract(objOrMat);

if isempty(targetTau) || ~isscalar(targetTau) || ~isfinite(targetTau) || targetTau <= 0
    error('escalation_advisor_20261022:target', ...
        'A positive scalar targetTau is required (mat input may supply it).')
end

W = iGet(opts, 'cellWeights', []);
if isempty(W)
    W = ones(size(slackMask));
else
    W = double(W);
    if ~isequal(size(W), size(slackMask))
        error('escalation_advisor_20261022:weights', ...
            'cellWeights size %s must match the grid %s.', ...
            mat2str(size(W)), mat2str(size(slackMask)))
    end
end

% --- Cell partition (the blocker fix lives here) ------------------------
% slack cells: solve clamped s = 0, Ftendon = 0, F_ins = 0 -> Torque_ins
%   identically 0 (construction, not physics: EXCLUDED from the torque
%   test, mapped to the RoM/slack side).
% infeasible cells: required contraction exceeds KMAX*Rest (the realized
%   travel deficit) -> RoM side, torque NaN there.
% taut feasible cells: the ONLY cells the torque-short test may read.
tauGrid = abs(TorqueIns);
nCells = numel(slackMask);
slackMask = slackMask & (W > 0);          % masked corners count nowhere
infeasMask = infeasMask & (W > 0);
tautMask = ~slackMask & ~infeasMask & isfinite(tauGrid) & (W > 0);

nSlack = nnz(slackMask);
nInfeas = nnz(infeasMask);
nTaut = nnz(tautMask);

% Raw (slack-gamed) margin over ALL finite cells -- reported, never used
% for the verdict. This is the number a slack-blind test would read.
validAll = tauGrid(~isnan(tauGrid));
if isempty(validAll)
    rawMargin = NaN;
else
    rawMargin = min(validAll) / targetTau - 1;
end

% Torque-short test: TAUT, FEASIBLE cells only.
if nTaut > 0
    tautTau = tauGrid(tautMask);
    [minTautTau, ix] = min(tautTau);
    tautMargin = minTautTau / targetTau - 1;
    idxCell = find(tautMask);
    minTautCell = iCellSubs(idxCell(ix), size(slackMask));
    torqueShort = minTautTau < targetTau;
    torqueEvaluable = true;
else
    minTautTau = NaN;
    tautMargin = NaN;
    minTautCell = [];
    torqueShort = false;
    torqueEvaluable = false;
end

% RoM/slack-side test: slack cells (over-fed travel / dead tendon) and
% infeasible cells (contraction demand beyond KMAX*Rest) are both
% travel-domain shortfalls; both step the tackle side (travel gain = G;
% higher G lowers the required contraction (DeltaL + delta_t)/G, which
% both un-slacks an over-fed span and shrinks the demand past KMAX).
romShort = (nSlack > 0) || (nInfeas > 0);

if nSlack > 0 && nInfeas > 0
    worstSlackSide = 'both';
elseif nSlack > 0
    worstSlackSide = 'slack';
elseif nInfeas > 0
    worstSlackSide = 'infeasible';
else
    worstSlackSide = '';
end

% --- Verdict and ladder step -------------------------------------------
reasons = {};
if torqueShort
    reasons{end + 1} = sprintf( ...
        ['torque-short: min taut-cell |Torque_ins| = %.4g N*m < target ' ...
         '%.4g N*m (taut margin %+.4f)'], minTautTau, targetTau, tautMargin);
else
    reasons{end + 1} = sprintf( ...
        'torque OK over %d taut feasible cell(s) (min %.4g vs target %.4g)', ...
        nTaut, minTautTau, targetTau);
end
if ~torqueEvaluable
    reasons{end + 1} = ['torque test NOT evaluable (no taut feasible ' ...
        'cells): cannot justify +BPA on a slack design'];
end
if nSlack > 0
    reasons{end + 1} = sprintf( ...
        ['RoM/slack side: %d of %d cells slack (Torque_ins = 0 there BY ' ...
         'CONSTRUCTION, excluded from the torque test)'], nSlack, nCells);
end
if nInfeas > 0
    reasons{end + 1} = sprintf( ...
        'RoM/slack side: %d of %d cells infeasible (travel demand > KMAX*Rest)', ...
        nInfeas, nCells);
end
if ~romShort && ~torqueShort
    reasons{end + 1} = 'no slack, no infeasible cells, taut torque meets target';
end

if torqueShort && romShort
    verdict = 'both';
elseif torqueShort
    verdict = 'torque_short';
elseif romShort
    verdict = 'rom_slack_short';
else
    verdict = 'converged';
end

% Ladder steps, capped at the directive grid {1,2,3} x {1,2,3}.
stepNBPA = 0; nextNBPA = NaN;
stepG = 0; nextG = NaN;
if torqueShort
    if curNBPA < capNBPA
        stepNBPA = 1;
        nextNBPA = curNBPA + 1;
        reasons{end + 1} = sprintf( ...
            'STEP: nBPA %d -> %d (insertion force gain = nBPA/G)', ...
            curNBPA, nextNBPA);
    else
        reasons{end + 1} = sprintf( ...
            'torque-short but nBPA = %d already at the directive cap %d', ...
            curNBPA, capNBPA);
    end
end
if romShort
    if curG < capG
        stepG = 1;
        nextG = curG + 1;
        reasons{end + 1} = sprintf( ...
            'STEP: G %d -> %d (travel gain = G shrinks required contraction)', ...
            curG, nextG);
    else
        reasons{end + 1} = sprintf( ...
            'RoM/slack-short but G = %d already at the directive cap %d', ...
            curG, capG);
    end
end
converged = strcmp(verdict, 'converged') || (stepNBPA == 0 && stepG == 0);

advice = struct( ...
    'verdict', verdict, ...
    'stepNBPA', stepNBPA, ...
    'stepG', stepG, ...
    'nextNBPA', nextNBPA, ...
    'nextG', nextG, ...
    'converged', converged, ...
    'reason', strjoin(reasons, '; '));

report = struct( ...
    'source', iSrcLabel(objOrMat, srcLabel), ...
    'targetTau', targetTau, ...
    'nCells', nCells, ...
    'nSlack', nSlack, ...
    'nInfeasible', nInfeas, ...
    'nTautFeasible', nTaut, ...
    'minTautTau', minTautTau, ...
    'tautMargin', tautMargin, ...
    'rawMargin', rawMargin, ...
    'minTautCell', minTautCell, ...
    'torqueTestEvaluable', torqueEvaluable, ...
    'worstSlackSide', worstSlackSide, ...
    'verdict', verdict, ...
    'curNBPA', curNBPA, 'curG', curG, ...
    'stepNBPA', stepNBPA, 'stepG', stepG);

if verbose
    fprintf(['--- escalation_advisor_20261022: %s | verdict = %s | ' ...
        'steps (nBPA,G) = (+%d,+%d) ---\n'], report.source, verdict, ...
        stepNBPA, stepG);
    for r = 1:numel(reasons)
        fprintf('    %s\n', reasons{r});
    end
    fprintf('    [raw slack-gamed margin %+.4f vs taut margin %s]\n', ...
        rawMargin, iNum(tautMargin));
end

% Optional CSV append (open/append/close).
if isfield(opts, 'csvFile') && ~isempty(opts.csvFile)
    fid = fopen(opts.csvFile, 'a');
    if fid > 0
        if ftell(fid) == 0
            fprintf(fid, ['source,targetTau,nCells,nSlack,nInfeasible,' ...
                'nTautFeasible,minTautTau,tautMargin,rawMargin,verdict,' ...
                'curNBPA,curG,stepNBPA,stepG\n']);
        end
        fprintf(fid, ['%s,%.6g,%d,%d,%d,%d,%.6g,%s,%s,%s,%d,%d,%d,%d\n'], ...
            report.source, targetTau, nCells, nSlack, nInfeas, nTaut, ...
            minTautTau, iNum(tautMargin), iNum(rawMargin), verdict, ...
            curNBPA, curG, stepNBPA, stepG);
        fclose(fid);
    end
end
end

% -------------------------------------------------------------------------
function [TorqueIns, slackMask, infeasMask] = iExtract(obj)
if isstruct(obj)
    % Synthetic/struct entry: grids supplied directly.
    req = {'TorqueIns', 'PulleySlack', 'PulleyInfeasible'};
    for k = 1:numel(req)
        if ~isfield(obj, req{k})
            fn = fieldnames(obj);
            fn = reshape(fn, 1, []);
            error('escalation_advisor_20261022:structInput', ...
                'Struct input needs fields %s (got %s).', ...
                strjoin(req, ', '), strjoin(fn, ', '))
        end
    end
    TorqueIns = obj.TorqueIns;
    slackMask = logical(obj.PulleySlack);
    infeasMask = logical(obj.PulleyInfeasible);
    if ~isequal(size(TorqueIns), size(slackMask)) || ...
            ~isequal(size(TorqueIns), size(infeasMask))
        error('escalation_advisor_20261022:structSize', ...
            'TorqueIns/Slack/Infeasible grids must share a size.')
    end
    return
end
cls = class(obj);
switch cls
    case 'BiPam_pulley'
        aK = obj.ActiveCrossing;
        slackMask = logical(obj.PulleySlack);
        infeasMask = logical(obj.PulleyInfeasible);
        % Torque_ins is N1 x 3 x N2 x 2 (per crossing); the scalar
        % torque magnitude about the active crossing's joint axis = the
        % norm over the 3-vector (the Opt_run_BiPulley margin uses
        % min over all finite |components|; the magnitude is the
        % configuration-diagnosis equivalent and is sign-safe).
        T3 = squeeze(obj.Torque_ins(:, :, :, aK));   % N1 x 3 x N2
        TorqueIns = reshape(vecnorm(T3, 2, 2), size(slackMask));
    case 'MonoPam_pulley'
        TorqueIns = obj.Torque_ins(:, 3);            % flexor z-torque
        slackMask = logical(obj.PulleySlack(:));
        infeasMask = logical(obj.PulleyInfeasible(:));
    otherwise
        error('escalation_advisor_20261022:input', ...
            'Unsupported input class %s (BiPam_pulley | MonoPam_pulley | struct | mat path).', cls)
end
if ~isequal(size(TorqueIns), size(slackMask)) || ...
        ~isequal(size(TorqueIns), size(infeasMask))
    error('escalation_advisor_20261022:size', ...
        'Torque/slack/infeasible grid size mismatch in %s input.', cls)
end
end

function v = iGet(opts, name, def)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = def;
end
end

function s = iNum(v)
if isnan(v)
    s = 'nan';
else
    s = sprintf('%.6g', v);
end
end

function s = iSrcLabel(obj, srcLabel)
if ~isempty(srcLabel)
    [~, f, e] = fileparts(srcLabel);
    s = [f, e];
else
    s = class(obj);
end
end

function c = iCellSubs(linIdx, szd)
if numel(szd) == 2
    [i, j] = ind2sub(szd, linIdx);
    c = [i, j];
else
    c = linIdx;
end
end
