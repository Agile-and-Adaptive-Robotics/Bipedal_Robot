%% smoke_Xi1gtXi2_20260912.m
% Pre-launch checks for the constrained (Xi1 > Xi2) minimizeFlxPin10mm rerun:
%  (1) checkcode the edited driver (gamultiobj call now carries nonlcon2)
%  (2) tiny constrained gamultiobj: call signature + constraint enforcement
%  (3) Xi1-vs-Xi2 violation stats on the unconstrained offT4 front

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(fileparts(here)));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));

%% (1) checkcode driver
msgs = checkcode('minimizeFlxPin10mm.m');
if isempty(msgs)
    fprintf('(1) checkcode: minimizeFlxPin10mm.m clean\n');
else
    fprintf('(1) checkcode: %d messages (first 5):\n', numel(msgs));
    for k = 1:min(5, numel(msgs))
        fprintf('    line %d: %s\n', msgs(k).line, msgs(k).message);
    end
end

%% (2) tiny constrained gamultiobj (serial, ~10 s)
opts = optimoptions('gamultiobj', 'UseParallel', false, 'Display', 'off', ...
    'PopulationSize', 12, 'MaxGenerations', 4);
[x, ~] = gamultiobj(@(X) sum(X, 2), 3, [], [], [], [], ...
    [0 log10(5e3) log10(5e3)], [2 log10(5e7) log10(5e7)], ...
    @(x) nlc(x), opts);
ok2 = all(x(:,2) >= x(:,3) - 1e-9);
fprintf('(2) constrained gamultiobj call form works; all %d front rows satisfy logXi1 >= logXi2: %s\n', ...
    size(x,1), string(ok2));

%% (3) violation stats on the unconstrained offT4 front
S = load('minimizeFlxPin10_results_20260911_2brkt_2trans_offT4.mat');
fr = S.filtered_results;   %[ind, hold(2), Xi0_m, Xi1, Xi2, train(3), val(3), dist]
xi1v = fr(:,5); xi2v = fr(:,6);
f1 = fr(:,2) == 1 & fr(:,3) == 2;   %fold-1 rows (Ben's split)
fprintf('(3) unconstrained offT4 front: %d/%d rows have Xi1 <= Xi2 (%.0f%%); fold-1 rows: %d/%d (%.0f%%)\n', ...
    sum(xi1v <= xi2v), numel(xi1v), 100*mean(xi1v <= xi2v), ...
    sum(xi1v(f1) <= xi2v(f1)), sum(f1), 100*mean(xi1v(f1) <= xi2v(f1)));
fprintf('    ratio Xi1/Xi2 across pooled front: min %.2f, median %.2f, max %.2f\n', ...
    min(xi1v./xi2v), median(xi1v./xi2v), max(xi1v./xi2v));

function [c, ceq] = nlc(x)
c = x(3) - x(2);   %same form as the driver's nonlcon2: enforces Xi2 < Xi1
ceq = [];
end
