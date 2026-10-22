%DEMO_KINEMATICS_GUARDS_20261022 Short live demonstration of the
%   2026-10-22 kinematics guards (plan Items 5 and 6). Runs both gate test
%   scripts in this folder and prints the envelope's control points and
%   untrusted banner (the Item 7 python audit has its own gate command).
%   Deterministic: no parpool, no optimizer calls, no MonoPam/BiPam
%   objects, no ModernRobotics path needed. Run:
%     matlab -batch "run('D:/GitHub/Bipedal_Robot/tmp/paper_bridge/opt_improve_20261022/demo_kinematics_guards_20261022.m')"

fprintf('==== kinematics guards demo (plan Items 5+6), %s ====\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

% Show the shipped envelope once, verbose (control points + banner)
fprintf('\n[Item 5] hip adduction/abduction coupling envelope, as shipped:\n');
[~, info] = hipAdductionEnvelope_20261022(0, 0, struct('verbose', true));
fprintf('trusted = %d (pending rows: %d of %d)\n\n', info.trusted, ...
    numel(info.pendingRows), height(info.table));

% Gate tests
test_hipAdductionEnvelope_20261022;
test_softCornerWeights_20261022;

fprintf('==== kinematics guards demo COMPLETE ====\n');
