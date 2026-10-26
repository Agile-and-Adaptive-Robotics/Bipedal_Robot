% Debug the demo-run failures from part A: pull the full error report.
root = 'D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape';
cd(root); addpath(root); addpath(fullfile(root,'demos'));
try
    run(fullfile(root,'demos','sns_run_demo.m'));
    disp('DEMO RUN OK');
catch e
    disp(getReport(e, 'extended', 'hyperlinks', 'off'));
end
