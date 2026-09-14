% patch_FlxPin10mats_allBPA_20260912.m
% Add lowercase allBPA + numBPA to every minimizeFlxPin10_results_202609*.mat
% (root + Dig_out archives). Non-destructive: save -append only adds/overwrites
% the two names; everything else (incl. ALLBPA) is left untouched.
% Priority: existing lowercase allBPA > ALLBPA > default 1:5.
% Ben, 2026-09-12: "modify the save results ... so that numBPA and allBPA appear".

here = fileparts(mfilename('fullpath'));
root = fileparts(here);   %2022_02_Festo
files = [ ...
    dir(fullfile(root, 'minimizeFlxPin10_results_202609*.mat')); ...
    dir(fullfile(here, 'old_bounds',   'minimizeFlxPin10_results_202609*.mat')); ...
    dir(fullfile(here, 'old_T3_results','minimizeFlxPin10_results_202609*.mat'))];

for iFile = 1:numel(files)
    fPath = fullfile(files(iFile).folder, files(iFile).name);
    S = load(fPath);
    hadLower = isfield(S, 'allBPA');
    hadNum   = isfield(S, 'numBPA');
    if hadLower
        allBPA = S.allBPA;
    elseif isfield(S, 'ALLBPA')
        allBPA = S.ALLBPA;
    else
        allBPA = 1:5;   %legacy vintages predate both names
    end
    allBPA = double(allBPA(:))';      %force a 1xN row of doubles
    numBPA = numel(allBPA);
    save(fPath, '-append', 'allBPA', 'numBPA');
    fprintf('%-60s allBPA(%s)=[%s] numBPA=%d  [had allBPA:%d numBPA:%d]\n', ...
        strrep(files(iFile).name, 'minimizeFlxPin10_results_', ''), ...
        ternary(hadLower, 'lower', ternary(isfield(S,'ALLBPA'), 'ALLBPA', 'default')), ...
        num2str(allBPA), numBPA, hadLower, hadNum);
end
fprintf('Done: %d mats patched.\n', numel(files));

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
