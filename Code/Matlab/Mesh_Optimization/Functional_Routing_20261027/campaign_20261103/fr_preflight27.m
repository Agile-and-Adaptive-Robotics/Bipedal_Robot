function ok = fr_preflight27()
%FR_PREFLIGHT27 Launch-safety pre-flight for run_0 (LAUNCH_DECISION.md).
%
% Builds ALL 27 contract-driven cases (fr_cases27), prints each case's
% structural summary, and evaluates the SEED design on BOTH branches
% (rigid + compliant Xi0-Xi3). NO solves here - this check exists to
% catch build/eval failures BEFORE the campaign launches. A group that
% fails gets QUEUED with its exact reason.
%
% Exit: ok = true iff every group built and seed-evaluated finitely on
% both branches. Diary: campaign output/preflight27.log.

C = fr_setup();
logFile = fullfile(C.outputDir, 'preflight27.log');
if exist(logFile, 'file')
    delete(logFile);
end
diary(logFile);
cleanupObj = onCleanup(@() diary('off'));

fprintf('=== preflight27 | %s | machine %s ===\n', datetime('now'), ...
    getenv('COMPUTERNAME'));
tAll = tic;
kin = fr_kinematics(C.osimStock);

failed = {};
ok = true;
try
    CASES = fr_cases27(C, kin);
catch err
    fprintf('CASE BUILD FAILED: %s | %s\n', err.identifier, err.message);
    ok = false;
    diary('off');
    error('fr_preflight27:build', 'fr_cases27 build failed: %s', ...
        err.message);
end
fprintf('cases built: %d | wall %.1f s\n', numel(CASES), toc(tAll));

evalWall = [];
for c = 1:numel(CASES)
    cs = CASES(c);
    tG = tic;
    fprintf(['\n[%2d] %-9s %-14s nRows %d cross %s nMov %d nBPA %d ' ...
        'dia %s tackle %d passes %d | target peak %.1f N*m (contract ' ...
        'mifN %.0f, %+.1f%%)\n'], cs.groupMapId, cs.group, cs.region, ...
        cs.nRows, mat2str(cs.cross), numel(cs.movable), ...
        cs.seed.nBPA, cs.diaStr, cs.hasTackle, numel(cs.passes), ...
        cs.target.peakAbs, cs.target.contractMifN, ...
        cs.target.peakVsContractPct);
    fprintf('     note: %s\n', strrep(cs.note, newline, ' '));
    bodySeq = strjoin(cellfun(@(p) p.body, cs.rows, ...
        'UniformOutput', false), '>');
    fprintf('     rows: %s\n', bodySeq);
    try
        oq = struct('quiet', true);
        [fR, repR] = fr_objective(C, cs, kin, cs.xSeed, 'rigid', oq);
        evalWall(end + 1) = repR.wallSeconds;
        [fC, repC] = fr_objective(C, cs, kin, cs.xSeed, 'compliant', oq);
        evalWall(end + 1) = repC.wallSeconds;
        finiteOk = isfinite(fR) && isfinite(fC) && ...
            all(isfinite(repR.tauPrimary)) && all(isfinite(repC.tauPrimary));
        fprintf(['     SEED rigid f %+.4f (viol %.2e) | compliant f ' ...
            '%+.4f (viol %.2e) | meet %d/%d poses (compliant) | ' ...
            'eval %.2f/%.2f s | %s\n'], fR, repR.violation, fC, ...
            repC.violation, sum(repC.meetsAtAllPoses), ...
            numel(repC.meetsAtAllPoses), repR.wallSeconds, ...
            repC.wallSeconds, fr_tf27(finiteOk));
        if isfield(cs, 'seedSpan') && ~strcmp(cs.seedSpan.kind, 'info')
            fprintf('     span %s: required %.4f m (allowance %.4f m)\n', ...
                cs.seedSpan.kind, cs.seedSpan.spanMax_m, ...
                cs.seedSpan.allowance_m);
        end
        if ~finiteOk
            ok = false;
            failed{end + 1} = sprintf('%s: nonfinite seed evaluation', ...
                cs.group);
        end
        if any(~isfinite([fR; fC]))
            ok = false;
            failed{end + 1} = sprintf('%s: NaN/Inf objective', cs.group);
        end
    catch err
        ok = false;
        failed{end + 1} = sprintf('%s: %s | %s', cs.group, ...
            err.identifier, err.message);
        fprintf('     SEED EVAL FAILED: %s | %s\n', err.identifier, ...
            err.message);
    end
end

fprintf('\n=== preflight27 wall %.1f min ===\n', toc(tAll) / 60);
evalWall = evalWall(isfinite(evalWall) & evalWall > 0);
if ~isempty(evalWall)
    fprintf(['seed-eval seconds: n %d | median %.3f | mean %.3f | ' ...
        'max %.3f\n'], numel(evalWall), median(evalWall), ...
        mean(evalWall), max(evalWall));
    nEvals = 27 * 2 * (C.solver.MaxFunctionEvaluations + ...
        C.refinementRounds * C.solver.refineMaxFunctionEvaluations);
    fprintf(['extrapolated full campaign at the seed-eval mean ' ...
        '(40500 evals): %.2f h serial\n'], nEvals * mean(evalWall) / 3600);
end
if ~ok
    fprintf('PREFLIGHT FAILED groups:\n  %s\n', strjoin(failed, '\n  '));
    error('fr_preflight27:failed', ...
        'preflight failed for %d group(s) - see %s', numel(failed), ...
        logFile);
end
fprintf('PREFLIGHT PASSED all 27 groups. log: %s\n', logFile);
end

function s = fr_tf27(b)
if b
    s = 'PASS';
else
    s = 'FAIL';
end
end
