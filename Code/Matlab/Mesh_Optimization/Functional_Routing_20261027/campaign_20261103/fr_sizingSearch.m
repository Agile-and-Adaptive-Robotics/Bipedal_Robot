function fr_sizingSearch(budgetEvals)
%FR_SIZINGSIZEARCH Discrete diameter/count/transmission screen (Ben
% directive item 4, 2026-11-04).
%
% Ben's ruling: allow 10 mm and 20 mm BPAs; the half-Fmax heuristic is the
% INITIAL SCREEN for choosing 10 mm, and the contract's preliminary 20 mm
% counts are search centers, NEVER mandatory. This search enumerates
% feasible diameter / nBPA / tackle-family / gain combinations per group
% and picks the winner by short-budget rigid-branch screen; the corrected
% campaign (run_1) then re-optimizes each winner config at FULL budget.
%
% Config grid per group:
%   nBPA    in {max(1, n-1), n, n+1}          (n = contract count, center)
%   gain    in {1, 2} for tackled groups      (1 = direct drive)
%   family  in {'L'} plus {'M'} for tackled   (catalog capacities apply)
%   dia     base from the 10 mm screen of record; the ALTERNATE diameter is
%           added only if its raw force-capacity screen is not hopeless:
%           peak|tauTarget| <= 0.9 * nBPA * Fmax_dia(at L0 seed)
%
% Screen: rigid branch, budget = C.solver.gateMaxFunctionEvaluations by
% default, seed start. Winner = best fBest among FEASIBLE configs
% (per-constraint tolerances, fr_objective v2); if none is feasible the
% best candidate is recorded with its top violator and run_1 handles the
% group explicitly (never silently dropped).
%
% Run headless: matlab -batch "addpath('<camp>'); fr_sizingSearch"
% Writes sizing_winners.json + output/sizing_search.log. SERIAL by design
% (safe beside other MATLAB work; no parpool contention with the
% campaign/agents). An optional numeric first argument overrides the
% per-config screen budget.

if nargin < 1
    budgetEvals = [];
end
budget = C0budget();
if isnumeric(budgetEvals) && ~isempty(budgetEvals)
    budget = double(budgetEvals(1));
end

logFile = fullfile('output', 'sizing_search.log');
if ~exist('output', 'dir'); mkdir('output'); end
fid = fopen(logFile, 'w');
cleanupFid = onCleanup(@() fclose(fid));
fprintf(fid, 'sizing search | %s | budget %d evals/config\n', ...
    datetime('now'), budget);

C = fr_setup();
if isfield(C, 'kinematics'); kin = C.kinematics; else; kin = fr_kinematics(C.osimStock); end
CASES = fr_cases27(C, kin);   % contract seeds; overrides applied per config
Fmax = FmaxByDia(C);

winners = struct('groups', {{}});
nDone = 0;
for c = 1:numel(CASES)
    cs = CASES(c);
    peakTau = max(abs(cs.target.tau));
    L0probe = cs.L0Seed;
    nBase = double(cs.seed.nBPA);
    diaBase = cs.diaStr;
    diaAlt = '20'; if strcmp(diaBase, '20'); diaAlt = '10'; end

    % ---- candidate configs ------------------------------------------------
    cfgs = cell(1, 32);   % preallocated, trimmed (max 2 dia x 3 n x 2 g x 2 fam)
    nCfg = 0;
    ns = unique(max(1, nBase + [-1, 0, 1]));
    diaList = {diaBase};
    cap = 0.9 * max(ns) * Fmax.(sprintf('d%s', diaBase))(L0probe);
    if peakTau <= cap
        diaList{2} = diaAlt;
    end
    for d = 1:numel(diaList)
        for n = ns
            if cs.hasTackle
                % a gain-1 tackle is mechanically contradictory (direct
                % drive wearing tackle hardware) and the class cannot
                % represent it - the sweep covers G=2 with both families;
                % direct-drive feasibility is the no-tackle groups' role
                for fam = {'L', 'M'}
                    nCfg = nCfg + 1;
                    cfgs{nCfg} = struct('diaStr', diaList{d}, ...
                        'nBPA', n, 'family', fam{1}, 'gain', 2);
                end
            else
                nCfg = nCfg + 1;
                cfgs{nCfg} = struct('diaStr', diaList{d}, ...
                    'nBPA', n, 'family', 'none', 'gain', 1);
            end
        end
    end
    cfgs = cfgs(1:nCfg);

    fprintf(fid, '\n=== %s (id %d) | base dia %s n %d tackle %d | peak|tauT| %.1f Nm | %d configs ===\n', ...
        cs.name, cs.groupMapId, diaBase, nBase, cs.hasTackle, peakTau, ...
        numel(cfgs));

    best = struct('fBest', inf, 'feasible', false, 'cfg', [], ...
        'topViolator', '', 'rep', []);
    for k = 1:numel(cfgs)
        ov = struct('diaStr', cfgs{k}.diaStr, 'nBPA', cfgs{k}.nBPA, ...
            'family', cfgs{k}.family, 'gain', cfgs{k}.gain);
        om = containers.Map('KeyType', 'char', 'ValueType', 'any');
        om(cs.name) = ov;
        csK = fr_cases27(C, kin, om);
        % select the INTENDED case (fr_cases27 returns all 27; csK(1) is
        % the contract's first group, NOT this group - gate caught the
        % mis-selection when every "winner" reproduced ercspn's score)
        selK = find(strcmp({csK.name}, cs.name), 1);
        assert(~isempty(selK), 'fr_sizingSearch:case', ...
            'override build lost group %s', cs.name);
        so = struct('quiet', true, 'budget', budget, 'envelopes', {{}});
        try
            [~, rep, info] = fr_solveCase(C, csK(selK), kin, 'rigid', so);
            score = info.fBest;
            fprintf(fid, '  cfg %2d: dia %s n %d fam %-2s g %d -> f %10.3f violEff %.2e feasible %d (%s)\n', ...
                k, cfgs{k}.diaStr, cfgs{k}.nBPA, cfgs{k}.family, ...
                cfgs{k}.gain, score, rep.violation, rep.feasible, ...
                rep.topViolator);
            better = (rep.feasible && ~best.feasible) ...
                || (rep.feasible == best.feasible && score < best.fBest);
            if better
                best = struct('fBest', score, 'feasible', rep.feasible, ...
                    'cfg', cfgs{k}, 'topViolator', rep.topViolator, ...
                    'rep', rep);
            end
        catch err
            st = getReport(err, 'extended', 'hyperlinks', 'off');
            stTail = st(max(1, numel(st) - 600):end);
            fprintf(fid, '  cfg %2d errored: %s | %s\n%s\n', k, ...
                err.identifier, err.message, stTail);
        end
    end

    if isempty(best.cfg) || ~isstruct(best.cfg)
        fprintf(fid, '  %s: NO config survived (all errored) - NOT-DONE for run_1\n', cs.name);
        winners.groups{end+1} = struct('name', cs.name, 'mapId', ...
            cs.groupMapId, 'status', 'no-config-survived');
       
    elseif best.feasible
        fprintf(fid, '  WINNER %s: dia %s nBPA %d fam %s gain %d, f %.3f (feasible at screen budget)\n', ...
            cs.name, best.cfg.diaStr, best.cfg.nBPA, best.cfg.family, ...
            best.cfg.gain, best.fBest);
        winners.groups{end+1} = struct('name', cs.name, 'mapId', ...
            cs.groupMapId, 'status', 'feasible', 'diaStr', ...
            best.cfg.diaStr, 'nBPA', best.cfg.nBPA, 'family', ...
            best.cfg.family, 'gain', best.cfg.gain, 'screenF', ...
            best.fBest);
    else
        fprintf(fid, '  BEST-INFEASIBLE %s: dia %s nBPA %d fam %s gain %d, f %.3f | %s (full budget in run_1)\n', ...
            cs.name, best.cfg.diaStr, best.cfg.nBPA, best.cfg.family, ...
            best.cfg.gain, best.fBest, best.topViolator);
        winners.groups{end+1} = struct('name', cs.name, 'mapId', ...
            cs.groupMapId, 'status', 'best-infeasible', 'diaStr', ...
            best.cfg.diaStr, 'nBPA', best.cfg.nBPA, 'family', ...
            best.cfg.family, 'gain', best.cfg.gain, 'screenF', ...
            best.fBest, 'topViolator', best.topViolator);
    end
    nDone = nDone + 1;
    fprintf('sizing: %d/%d groups screened (%s)\n', nDone, ...
        numel(CASES), cs.name);
end

outPath = fullfile('sizing_winners.json');
fid2 = fopen(outPath, 'w');
fwrite(fid2, jsonencode(winners));
fclose(fid2);
fprintf(fid, '\nwinners written: %s\n', outPath);
fprintf('sizing search complete: %s\n', outPath);
end

% ==========================================================================
function b = C0budget()
C = fr_setup();
b = C.solver.gateMaxFunctionEvaluations;
end

function F = FmaxByDia(C)
% zero-contraction capacity at the 620 kPa pressure of record, evaluated
% at the probe length by the repo's own force law (maxBPAforce.m:
% maxBPAforce(restingLength_m, diameterStr, pressure_kPa))
F.d10 = @(L) maxBPAforce(L, '10', C.pressure);
F.d20 = @(L) maxBPAforce(L, '20', C.pressure);
end
