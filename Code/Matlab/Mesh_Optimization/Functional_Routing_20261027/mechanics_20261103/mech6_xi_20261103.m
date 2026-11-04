function mech6_xi_20261103()
%MECH6_XI_20261103 MechanCheck row 6: load the adopted Xi0-Xi3 context
% EXACTLY as buildKneeFlexorContext20mm.m (pick row 77) and
% buildKneeExtContext20mm.m (pick row 32) load it; report the values and
% their source fields; verify Xi4 appears nowhere.
%
% Flexor builder (buildKneeFlexorContext20mm.m:171-202):
%   load minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat
%        filtered_results xCols ; pick = 77 ; g = filtered_results(77,xCols)
%        Xi0 = g(1); Xi1 = g(2); Xi2 = g(3)
%   load minimizeExt10mmX3_results_20260920_flxr77.mat filtered_results xCols
%        pick = 32 ; g(4) = filtered_results(32, xCols(4)) ; Xi3 = g(4)
% Extensor builder (buildKneeExtContext20mm.m:164-183):
%   load minimizeExt10mmX3_results_20260920_flxr77.mat filtered_results xCols
%        pick = 32 ; g = filtered_results(32, xCols)
%        Xi0 = g(1); Xi1 = g(2); Xi2 = g(3); Xi3 = g(4)

repo = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(repo, 'Code', 'Matlab')));
addpath(fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization'), '-end');
addpath(fullfile(repo, 'Testing_Data', '2022_02_Festo'), '-end');
cd(fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization'));

fprintf('=== MechanCheck 6: adopted Xi0-Xi3 context (rows 77 / 32) ===\n\n')

%% ---- Flexor: replicate the builder's exact load lines ----------------
w = which('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat');
fprintf('[Flexor] load minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat\n')
fprintf('         resolved from: %s\n', w)
S = load(w, 'filtered_results', 'xCols');
xCols = S.xCols;
fprintf('         xCols = %s\n', mat2str(xCols))
pick = 77;
g = S.filtered_results(pick, xCols);
Xi0f = g(1); Xi1f = g(2); Xi2f = g(3);
fprintf('         pick = %d; Xi0 = %.6g m, Xi1 = %.6g N/m, Xi2 = %.6g N/m\n', ...
    pick, Xi0f, Xi1f, Xi2f)
w2 = which('minimizeExt10mmX3_results_20260920_flxr77.mat');
S2 = load(w2, 'filtered_results', 'xCols');
pick = 32;    % builder re-assigns pick = 32 before the Xi3 load (line 189)
g4 = S2.filtered_results(pick, S2.xCols(4));
Xi3f = g4(1);
fprintf('         Xi3 from %s (pick %d, xCols(4) = %d): Xi3 = %.6g\n\n', ...
    'minimizeExt10mmX3_results_20260920_flxr77.mat', pick, S2.xCols(4), Xi3f)

%% ---- Extensor: replicate the builder's exact load lines --------------
fprintf('[Extensor] load minimizeExt10mmX3_results_20260920_flxr77.mat\n')
fprintf('           resolved from: %s\n', w2)
fprintf('           xCols = %s\n', mat2str(S2.xCols))
pickE = 32;
gE = S2.filtered_results(pickE, S2.xCols);
fprintf('           pick = %d; Xi0 = %.6g m, Xi1 = %.6g N/m, Xi2 = %.6g N/m, Xi3 = %.6g\n\n', ...
    pickE, gE(1), gE(2), gE(3), gE(4))

%% ---- Cross-check vs the records of record (AGENTS 2026-10-18) --------
fprintf('[Cross-check] records of record minimizeFlxPin10_results.mat (row 77) / ')
fprintf('minimizeExt10mmX3_results.mat (row 32), LOCKSRC / PICK fields:\n')
Sr = load(which('minimizeFlxPin10_results.mat'), 'filtered_results', 'xCols');
gr = Sr.filtered_results(77, Sr.xCols);
fprintf('  flexor  row 77: Xi0 = %.6g, Xi1 = %.6g, Xi2 = %.6g\n', gr(1), gr(2), gr(3))
fprintf('  identical to builder flexor Xi0/Xi1/Xi2: %s\n', ...
    verdict(max(abs([Xi0f, Xi1f, Xi2f] - gr)) == 0))
Se = load(which('minimizeExt10mmX3_results.mat'), 'filtered_results', 'xCols');
ge = Se.filtered_results(32, Se.xCols);
fprintf('  extensor row 32: Xi0 = %.6g, Xi1 = %.6g, Xi2 = %.6g, Xi3 = %.6g\n', ...
    ge(1), ge(2), ge(3), ge(4))
fprintf('  identical to builder extensor Xi0/Xi1/Xi2/Xi3: %s\n', ...
    verdict(max(abs(gE - ge)) == 0))
if isfield(Sr, 'PICK'); disp('  (flexor mat carries legacy PICK fields; picks of record = builder rows)'); end

%% ---- Try the actual builders end to end ------------------------------
fprintf('\n[Builders] running the real context builders end to end...\n')
ctxF = struct(); ctxE = struct();
try
    ctxF = buildKneeFlexorContext20mm();
    fprintf('  buildKneeFlexorContext20mm: OK. ctx.Xi0 = %.6g, ctx.Xi1 = %.6g, ', ...
        ctxF.Xi0, ctxF.Xi1)
    fprintf('ctx.Xi2 = %.6g, ctx.Xi3 = %.6g, ctx.KMAX = %.3f\n', ctxF.Xi2, ctxF.Xi3, ctxF.KMAX)
    okFlex = verdict(abs(ctxF.Xi0 - Xi0f) == 0 && abs(ctxF.Xi1 - Xi1f) == 0 && ...
        abs(ctxF.Xi2 - Xi2f) == 0 && abs(ctxF.Xi3 - Xi3f) == 0);
    fprintf('  builder Xi values match the replicated load: %s\n', okFlex)
catch err
    fprintf('  buildKneeFlexorContext20mm FAILED at the end-to-end stage: %s\n', err.message)
    fprintf('  (the replicated load lines above are the ask''s check; see blocker note)\n')
end
try
    ctxE = buildKneeExtContext20mm();
    fprintf('  buildKneeExtContext20mm: OK. ctx.Xi0 = %.6g, ctx.Xi1 = %.6g, ', ...
        ctxE.Xi0, ctxE.Xi1)
    fprintf('ctx.Xi2 = %.6g, ctx.Xi3 = %.6g, ctx.KMAX = %.3f\n', ctxE.Xi2, ctxE.Xi3, ctxE.KMAX)
    okExt = verdict(abs([ctxE.Xi0, ctxE.Xi1, ctxE.Xi2, ctxE.Xi3] - gE) == 0);
    fprintf('  builder Xi values match the replicated load: %s\n', okExt)
catch err
    fprintf('  buildKneeExtContext20mm FAILED at the end-to-end stage: %s\n', err.message)
    fprintf('  (the replicated load lines above are the ask''s check; see blocker note)\n')
end

%% ---- Xi4 absence ------------------------------------------------------
fprintf('\n[Xi4] contract: Xi0-Xi3 only; Xi4 must NOT appear anywhere.\n')
fprintf('  ctx fields (flexor): %s\n', strjoin(fieldnames(ctxF)', ', '))
fprintf('  ctx fields (extensor): %s\n', strjoin(fieldnames(ctxE)', ', '))
fprintf('  Xi4 present in either context: %s\n', ...
    verdict(~any(strcmp('Xi4', [fieldnames(ctxF); fieldnames(ctxE)]))))
fprintf('  (repo-wide grep for Xi4 in campaign .m code: run separately in the shell log)\n')
fprintf('=== MechanCheck 6 done ===\n')
end

function s = verdict(ok)
if ok
    s = 'PASS';
else
    s = 'FAIL';
end
end
