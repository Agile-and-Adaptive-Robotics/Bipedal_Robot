% rendercheck_FlxPin10mm_20260912.m
% Reproduce Ben's flow: load the full-workspace result mat, then run the
% driver's plotting sections (extracted verbatim into _plotsections_..._body.m),
% and export every figure to PNG for visual QA of the new panel letters + legend rule.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
cd(root);
addpath(root); addpath(here);
addpath(fullfile(root, '..', '..', 'Code', 'Matlab', 'Functions')); %ForceStrainForFit.mat lives here

matFile = fullfile(root, 'minimizeFlxPin10_results_20260912_2brkt_2trans_K2allX1_2folds.mat');
try
    S = load(matFile);   %full workspace incl. bpa, allBPA, numBPA
    fn = fieldnames(S);
    for k = 1:numel(fn), eval([fn{k} ' = S.(fn{k});']); end %#ok<SAGROW>
    fprintf('Full-workspace load OK (%d vars).\n', numel(fn));
catch ME
    fprintf('Full load failed (%s) -- loading bpa/allBPA/numBPA only.\n', ME.message);
    L = load(matFile, 'bpa', 'allBPA', 'numBPA');
    bpa = L.bpa; allBPA = L.allBPA; numBPA = L.numBPA;
end

run(fullfile(here, 'plotsectionsBody_FlxPin10mm_20260912.m'));

figs = struct('fig', {figTpre, figTpost, figL, figMA, figS}, ...
              'name', {'TorquePre','TorquePost','MuscleLength','MomentArm','RelStrain'});
for k = 1:numel(figs)
    outFile = fullfile(here, sprintf('plotcheck_FlxPin10mm_%s_20260912.png', figs(k).name));
    exportgraphics(figs(k).fig, outFile, 'Resolution', 110);
    fprintf('wrote %s\n', outFile);
end
close all force;
disp('RENDERCHECK DONE');
