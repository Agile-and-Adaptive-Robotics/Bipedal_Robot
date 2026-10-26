function run_desktop_capability_20261024(part)
% Desktop (EB475WS4, R2025a) capability pass over SNS_Simscape, 2026-10-24.
% Purpose: while the laptop (the only Simscape Multibody seat) is down,
% establish what of the base-Simulink SNS stack runs on THIS machine, per
% README_SNS_Simscape.md's documented run order and PASS numbers.
% Regeneration overwrites the laptop's R2025b-native SNS_Library.slx and
% demo .slx files with R2025a-native copies (documented path for R2025a
% machines); pre-pass backups live in dev\backup_pre_20261024\.
% Parts: 'a' env+library+tests+demos | 'b' deng+spinal+winner+mex | 'c' multibody load probe.

root = fileparts(fileparts(mfilename('fullpath')));   % dev\ -> SNS_Simscape
cd(root);
addpath(root); addpath(fullfile(root,'demos'));
logf = fullfile('logs', sprintf('desktop_capability_20261024_%s.log', part));
if exist(logf,'file'); delete(logf); end
diary(logf);
fprintf('=== SNS_Simscape desktop capability pass | host %s | %s | part %s ===\n', ...
        getenv('COMPUTERNAME'), version, part);
fprintf('cwd: %s\n', pwd);

switch part
case 'a'
    v = ver;
    fprintf('--- installed Simulink/Simscape products ---\n');
    for k = 1:numel(v)
        nm = v(k).Name;
        if contains(nm,'Simulink','IgnoreCase',true) || ...
           contains(nm,'Simscape','IgnoreCase',true)
            fprintf('%s | %s\n', nm, v(k).Version);
        end
    end
    step('sns_build_library',    fullfile(root,'sns_build_library.m'));
    step('sns_build_actuators',  fullfile(root,'sns_build_actuators.m'));
    step('sns_test_actuators',   fullfile(root,'sns_test_actuators.m'));
    step('sns_units_test_2n',    fullfile(root,'sns_units_test_2n.m'));
    step('sns_build_demo',       fullfile(root,'demos','sns_build_demo.m'));
    step('sns_run_demo',         fullfile(root,'demos','sns_run_demo.m'));
    step('sns_build_cpg_demo',   fullfile(root,'demos','sns_build_cpg_demo.m'));
    step('sns_run_cpg_demo',     fullfile(root,'demos','sns_run_cpg_demo.m'));
    step('sns_build_beer_demo',  fullfile(root,'demos','sns_build_beer_demo.m'));
    step('sns_run_beer_demo',    fullfile(root,'demos','sns_run_beer_demo.m'));
case 'b'
    step('sns_build_deng_cpg',      fullfile(root,'demos','sns_build_deng_cpg.m'));
    step('sns_run_deng_demo',       fullfile(root,'demos','sns_run_deng_demo.m'));
    step('sns_build_from_json',     fullfile(root,'sns_build_from_json.m'));
    step('sns_verify_from_json',    fullfile(root,'sns_verify_from_json.m'));
    step('sns_build_winner_slice',  fullfile(root,'demos','sns_build_winner_slice.m'));
    step('sns_run_winner_slice',    fullfile(root,'demos','sns_run_winner_slice.m'));
    fprintf('--- mex C compiler setup ---\n');
    try
        mex -setup C;
        fprintf('MEX_C_SETUP OK\n');
    catch e
        fprintf('MEX_C_SETUP FAIL: %s\n', e.message);
    end
case 'c'
    % Load-only probe of an imported Simscape Multibody model (no update,
    % no sim). Expected: license error (SimMechanics seats reserved).
    % Kept in its own process: README documents the crash-at-exit hazard
    % after a FAILED Simscape update.
    mfile = 'mdl_leg_rig_ba003_imported_R2025a.slx';
    try
        m = load_system(fullfile(root, mfile));
        fprintf('MULTIBODY_LOAD OK (load only, no update/sim): %s\n', m);
        close_system(m, 0);
    catch e
        fprintf('MULTIBODY_LOAD FAIL: %s\n', e.message);
    end
otherwise
    error('unknown part %s', part);
end
diary off;
fprintf('=== part %s done ===\n', part);
end

function step(name, scriptPath)
fprintf('--- STEP %s ---\n', name);
t0 = tic;
try
    run(scriptPath);
    fprintf('STEP %s OK (%.1f s)\n', name, toc(t0));
catch e
    fprintf('STEP %s FAIL (%.1f s): %s\n', name, toc(t0), e.message);
end
end
