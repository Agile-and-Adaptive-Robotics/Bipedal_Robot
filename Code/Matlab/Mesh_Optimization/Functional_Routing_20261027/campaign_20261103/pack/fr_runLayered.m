function result = fr_runLayered(C, kin, cases, opts)
%FR_RUNLAYERED Layered, iterative packaging driver (master item 5,
% Ben's ruling).
%
% PASS 1 (deep-outward): groups optimize DEEP-OUTWARD in the contract's
% layerOrder. Each finished layer's occupied envelope (sleeve + fittings
% + tackle + hardware, pose-dependent, fr_envelope) is RETAINED and
% enforced as a HARD interference obstacle inside every more-superficial
% group's constraint set (fr_objective constraint 7), so superficial
% muscles must route around what is already placed. The deepest layer
% optimizes against bone only.
%
% REFINEMENT ROUNDS (opts.rounds, default C.refinementRounds): any layer
% may re-solve while its neighbors' CURRENT envelopes constrain it, so an
% early deep choice cannot permanently block a superficial muscle. A
% round is ACCEPTED only if it does not degrade any retained layer's
% optimized objective beyond C.refineTolerance; otherwise that layer's
% previous design is kept (revert rule) and the non-degradation is
% recorded.
%
% result.placed(mapId): design, rep, envelope, objective-of-record.
if nargin < 4
    opts = struct();
end
if ~isfield(opts, 'rounds'); opts.rounds = C.refinementRounds; end
if ~isfield(opts, 'budget'); opts.budget = C.solver.MaxFunctionEvaluations; end
if ~isfield(opts, 'quiet'); opts.quiet = false; end

byId = containers.Map('KeyType', 'double', 'ValueType', 'any');
for c = 1:numel(cases)
    byId(cases(c).groupMapId) = c;
end
sched = fr_layers(C, cases);
result = struct();
result.placed = containers.Map('KeyType', 'double', 'ValueType', 'any');
result.log = cell(numel(sched.scheduledIds) * (1 + opts.rounds), 1);
logCount = 0;
envelopes = cell(1, numel(sched.scheduledIds));
envCount = 0;

% ---- PASS 1: deep-outward ------------------------------------------------
for s = 1:numel(sched.scheduledIds)
    id = sched.scheduledIds(s);
    cname = sched.caseByMapId(id);
    cs = cases(byId(id));
    so = opts; so.envelopes = envelopes;
    [best, rep, info] = fr_solveCase(C, cs, kin, 'rigid', so);
    rep.fRecord = info.fBest;
    env = fr_envelope(C, cs, best, kin, ...
        cs.passes{cs.envelopePassIdx}.qStructs, 'rigid');
    result.placed(id) = struct('design', best, 'rep', rep, ...
        'envelope', env, 'caseName', cname, 'info', info, ...
        'round', 0);
    envCount = envCount + 1; envelopes{envCount} = env;
    logCount = logCount + 1;
    result.log{logCount} = sprintf( ...
        'PASS1 %s (id %d, %s): f %.4f (viol %.2e)', cname, id, ...
        sched.region{s}, info.fBest, rep.violation);
end

% ---- REFINEMENT ROUNDS ----------------------------------------------------
for rd = 1:opts.rounds
    roundLog = cell(numel(sched.scheduledIds), 1);
    rlCount = 0;
    acceptedAll = true;
    for s = 1:numel(sched.scheduledIds)
        id = sched.scheduledIds(s);
        placed = result.placed(id);
        cs = cases(byId(id));
        % neighbors' envelopes = every OTHER placed layer's CURRENT envelope
        envOthers = {};
        fPrev = placed.rep.fRecord;
        placedIds = cellfun(@str2double, result.placed.keys());
        for id2 = placedIds
            if id2 == id; continue; end
            envOthers{end+1} = result.placed(id2).envelope; %#ok<AGROW>
        end
        so = opts; so.envelopes = envOthers;
        so.budget = C.solver.refineMaxFunctionEvaluations;
        [bestNew, repNew, infoNew] = fr_solveCase(C, cs, kin, ...
            'rigid', so);
        % accept only if not worse than the retained objective
        if infoNew.fBest <= fPrev + C.refineTolerance
            envNew = fr_envelope(C, cs, bestNew, kin, ...
                cs.passes{cs.envelopePassIdx}.qStructs, 'rigid');
            result.placed(id) = struct('design', bestNew, ...
                'rep', repNew, 'envelope', envNew, ...
                'caseName', placed.caseName, 'info', infoNew, ...
                'round', rd);
            rlCount = rlCount + 1;
            roundLog{rlCount} = sprintf( ...
                'REFINE r%d %s (id %d): f %.4f -> %.4f (accepted)', ...
                rd, placed.caseName, id, fPrev, infoNew.fBest);
        else
            acceptedAll = false;
            rlCount = rlCount + 1;
            roundLog{rlCount} = sprintf( ...
                ['REFINE r%d %s (id %d): f %.4f -> %.4f REJECTED ' ...
                '(degrades the retained objective > %.1e; kept r%d)'], ...
                rd, placed.caseName, id, fPrev, infoNew.fBest, ...
                C.refineTolerance, placed.round);
        end
    end
    logCount = logCount + 1;
    result.log{logCount} = strjoin(roundLog, ' | ');
    result.roundAccepted(rd) = acceptedAll;
end
end
