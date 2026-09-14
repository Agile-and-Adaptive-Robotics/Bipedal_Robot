%% Dig_a0shiftCheck_20260911.m
% How does the per-test RIGID baseline (a0) respond to the +5.3 deg encoder
% shift? Reuses the current evaluator's own rigid prediction
% (minimizeFlxPin(0,Inf,Inf)) and recomputes GoF at shifted experimental
% angles. Read-only: touches no file the offT4 CV run depends on.

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

[~, b0] = minimizeFlxPin(0, Inf, Inf);   %rigid baseline structs

% applied shift in the CURRENT evaluator build: +5.3 on test 4 only
applied = [0 0 0 5.3 0];
tests = [3 4];
shifts = [-10.6 -5.3 -2.65 0 2.65 5.3 10.6];
lbl = {'48cm','46cm','47cm','40cm-t','41cm'};

for j = tests
    [Ak_s, idx] = sort(b0(j).Ak);
    M_s = b0(j).M_p(idx, 3);
    F = griddedInterpolant(Ak_s, M_s);
    raw = b0(j).Aexp - applied(j);    %raw reported experimental angles
    fprintf('\nTest %d (%s) -- rigid baseline vs angle shift on the experimental angles:\n', j, lbl{j});
    fprintf('  shift(deg)   RMSE      FVU      MaxRes\n');
    for s = shifts
        Mopt = F(raw + s);
        [r, fv, mr] = Go_OfF(b0(j).Mexp, Mopt);
        fprintf('  %7.2f  %7.3f  %7.4f  %8.3f\n', s, r, fv, mr);
    end
end
fprintf('\n(self-check: shift = applied must reproduce the a0 table in the run log)\n');
