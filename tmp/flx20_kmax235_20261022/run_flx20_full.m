%% run_flx20_full.m
% Runs Ben's Knee_Flexor_Data_20mm.m (edited 2026-10-22: KMAX 0.235)
% from this tmp dir so its cwd-written outputs (GIF, CSV) stay out of the
% repo trees. The OpenSim_Bifem txts are copied here for the relative
% readmatrix calls. Figures are saved as PNGs for review.
% NOTE: the run script starts with `clear`, which wipes this wrapper's
% workspace -- re-derive everything AFTER run() returns or throws.
clc; clear; close all;
here = fileparts(mfilename('fullpath'));
cd(here);
try
    run('D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Results\Knee_Flexor_Data_20mm.m');
    runOK = true;
catch e
    runOK = false;
    errTxt = getReport(e, 'basic', 'hyperlinks', 'off');
end
here = fileparts(mfilename('fullpath'));
cd(here);
if runOK
    fprintf('\n[FULL RUN OK]\n');
else
    fprintf('\n[FULL RUN ERROR]\n%s\n', errTxt);
end
figs = orderdoc(findall(0, 'Type', 'figure'));
for i = 1:numel(figs)
    nm = get(figs(i), 'Name');
    if isempty(nm), nm = sprintf('fig%02d', i); end
    fn = matlab.lang.makeValidName(nm);
    saveas(figs(i), fullfile(here, sprintf('%02d_%s.png', i, fn)));
end
fprintf('[SAVED %d FIGURE PNGs to %s]\n', numel(figs), here);
