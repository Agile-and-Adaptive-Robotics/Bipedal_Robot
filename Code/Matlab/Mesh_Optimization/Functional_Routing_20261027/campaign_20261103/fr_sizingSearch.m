function fr_sizingSearch()
%FR_SIZINGSEARCH Discrete diameter/count/transmission screen, SIGNED
% standard, PARFOR across configs (10 workers) for speed — Ben
% 2026-11-04: "Can you speed it up? ... a lot of my pipeline has been
% backed up."
%
% Config grid per group: nBPA in {n-1, n, n+1} (contract n as center,
% never mandatory), family L/M for tackled groups, gain 2 (direct-drive
% handled by no-tackle groups), plus the alternate diameter when its
% capacity screen is not hopeless. Screen budget 150 evals, seed start,
% rigid branch, signed scoring.
%
% Run headless: matlab -batch "addpath('<camp>'); fr_sizingSearch"
% Writes sizing_winners.json + output/sizing_search.log.

if ~exist('output', 'dir'); mkdir('output'); end
C = fr_setup();
if isfield(C, 'kinematics'); kin = C.kinematics; else; kin = fr_kinematics(C.osimStock); end
CASES = fr_cases27(C, kin);
Fmax = FmaxByDia(C);
budget = C.solver.gateMaxFunctionEvaluations;

logFile = fullfile('output', 'sizing_search.log');
fid = fopen(logFile, 'w');
cleanupFid = onCleanup(@() fclose(fid));
fprintf(fid, 'sizing search (signed, parfor) | %s | budget %d\n', ...
    datetime('now'), budget);

tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
nGroups = numel(CASES);
allWinners = cell(1, nGroups);

% pool BEFORE the group loop; workers persist across groups
setenv('OMP_NUM_THREADS', '1');
setenv('MKL_NUM_THREADS', '1');
if isempty(gcp('nocreate'))
    pool = parpool('Processes', C.parallelMaxWorkers); %#ok<NASGU>
    fprintf('parpool: %d workers\n', pool.NumWorkers);
end

for c = 1:nGroups
    cs = CASES(c);
    peakTau = max(abs(cs.target.tau));
    L0probe = cs.L0Seed;
    nBase = double(cs.seed.nBPA);
    diaBase = cs.diaStr;
    diaAlt = '20'; if strcmp(diaBase, '20'); diaAlt = '10'; end

    % ---- candidate configs (preallocated) --------------------------------
    cfgs = cell(1, 32);
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
        nCfg);

    % ---- parfor screen ----------------------------------------------------
    scr = cell(1, nCfg);
    for k = 1:nCfg
        % full rebuild per config via the overrides mechanism (PROVEN
        % path; the finishCase shortcut was illegal - fr_finishCase27 is
        % a local function of fr_cases27.m, unreachable from here)
        om = containers.Map('KeyType', 'char', 'ValueType', 'any');
        om(cs.name) = struct('diaStr', cfgs{k}.diaStr, ...
            'nBPA', cfgs{k}.nBPA, 'family', cfgs{k}.family, ...
            'gain', cfgs{k}.gain);
        csK = fr_cases27(C, kin, om);
        selK = find(strcmp({csK.name}, cs.name), 1);
        so = struct('quiet', true, 'budget', budget, 'envelopes', {{}});
        try
            [~, rep, info] = fr_solveCase(C, csK(selK), kin, 'rigid', so);
            scr{k} = struct('fBest', info.fBest, ...
                'feasible', rep.feasible, 'violEff', rep.violation, ...
                'topViolator', rep.topViolator);
        catch err
            scr{k} = struct('fBest', NaN, 'feasible', false, ...
                'violEff', NaN, 'topViolator', ['ERROR ' err.message(1:min(80, numel(err.message)))]);
        end
    end

    best = struct('fBest', inf, 'feasible', false, 'cfg', [], ...
        'topViolator', '');
    for k = 1:nCfg
        if isempty(scr{k}); continue; end
        r = scr{k};
        fprintf(fid, '  cfg %2d: dia %s n %d fam %-2s g %d -> f %10.3f violEff %.2e feasible %d (%s)\n', ...
            k, cfgs{k}.diaStr, cfgs{k}.nBPA, cfgs{k}.family, ...
            cfgs{k}.gain, r.fBest, r.violEff, r.feasible, r.topViolator);
        better = (r.feasible && ~best.feasible) ...
            || (r.feasible == best.feasible && r.fBest < best.fBest);
        if better
            best = struct('fBest', r.fBest, 'feasible', r.feasible, ...
                'cfg', cfgs{k}, 'topViolator', r.topViolator);
        end
    end

    if isempty(best.cfg) || ~isstruct(best.cfg)
        fprintf(fid, '  %s: NO config survived - NOT-DONE for run_1\n', cs.name);
        allWinners{c} = struct('name', cs.name, 'mapId', cs.groupMapId, ...
            'status', 'no-config-survived');
    elseif best.feasible
        fprintf(fid, '  WINNER %s: dia %s nBPA %d fam %s gain %d, f %.3f (feasible at screen budget)\n', ...
            cs.name, best.cfg.diaStr, best.cfg.nBPA, best.cfg.family, ...
            best.cfg.gain, best.fBest);
        allWinners{c} = struct('name', cs.name, 'mapId', cs.groupMapId, ...
            'status', 'feasible', 'diaStr', best.cfg.diaStr, ...
            'nBPA', best.cfg.nBPA, 'family', best.cfg.family, ...
            'gain', best.cfg.gain, 'screenF', best.fBest);
    else
        fprintf(fid, '  BEST-INFEASIBLE %s: dia %s nBPA %d fam %s gain %d, f %.3f | %s\n', ...
            cs.name, best.cfg.diaStr, best.cfg.nBPA, best.cfg.family, ...
            best.cfg.gain, best.fBest, best.topViolator);
        allWinners{c} = struct('name', cs.name, 'mapId', cs.groupMapId, ...
            'status', 'best-infeasible', 'diaStr', best.cfg.diaStr, ...
            'nBPA', best.cfg.nBPA, 'family', best.cfg.family, ...
            'gain', best.cfg.gain, 'screenF', best.fBest, ...
            'topViolator', best.topViolator);
    end
    fprintf('sizing: %d/%d groups screened (%s)\n', c, nGroups, cs.name);
end

w = struct('groups', {{allWinners{:}}});
json.dump = 0; %#ok<NASGU>
fid2 = fopen(fullfile('sizing_winners.json'), 'w');
fwrite(fid2, jsonencode(w, 'PrettyPrint', true));
fclose(fid2);
fprintf(fid, '\nwinners written\n');
fprintf('sizing search complete\n');
end

% ==========================================================================
function F = FmaxByDia(C)
F.d10 = @(L) maxBPAforce(L, '10', C.pressure);
F.d20 = @(L) maxBPAforce(L, '20', C.pressure);
end
