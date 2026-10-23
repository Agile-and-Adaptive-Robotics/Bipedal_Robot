% Check Verbatim 20261023
% Regression check: the 4 routes the original biPulleySpecsFromMasters.m
% carries (semimem_r, rect_fem_r, med_gas_r, tfl_r) must come out of the
% dated extension biPulleySpecsFromMasters_20261023 IDENTICALLY (points,
% cross, frames, diameter, pulley config, fmax, and every kin field the
% original emits) when both builders read the SAME masters data (the
% campaign's read-only repo copies at tmp/old_code_copy_20261022).
% The original file itself is NOT modified; it is CALLED with mastersDir
% pointed at the repo copies (its default thumb-drive path may not exist
% on this machine).

function check_verbatim_20261023()

scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    root = parent;
end
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(scriptDir);

dataDir = fullfile(root, 'tmp', 'old_code_copy_20261022');
shared = {'semimem_r', 'rect_fem_r', 'med_gas_r', 'tfl_r'};

orig = biPulleySpecsFromMasters(shared, dataDir);
new = biPulleySpecsFromMasters_20261023(shared, dataDir);

nDiff = 0;
for k = 1:numel(shared)
    A = orig(k);  B = new(k);
    assert(strcmp(A.name, B.name));
    checks = { ...
        'points',        isequal(A.points, B.points); ...
        'cross',         isequal(A.cross, B.cross); ...
        'frames',        isequal(A.frames, B.frames); ...
        'diameter',      isequal(A.diameter, B.diameter); ...
        'pulley',        isequal(A.pulley, B.pulley); ...
        'bundleOffsets', isequal(A.bundleOffsets, B.bundleOffsets); ...
        'bundleMargin',  isequal(A.bundleMargin, B.bundleMargin); ...
        'footprintWidth', isequal(A.footprintWidth, B.footprintWidth); ...
        'fmax',          isequal(A.fmax, B.fmax); ...
        'fmaxSource',    isequal(A.fmaxSource, B.fmaxSource); ...
        'kin.pivots',    isequal(A.kin.pivots, B.kin.pivots); ...
        'kin.axes',      isequal(A.kin.axes, B.kin.axes); ...
        'kin.types',     isequal(A.kin.types, B.kin.types); ...
        'kin.thetaRangesDeg', isequal(A.kin.thetaRangesDeg, ...
                                      B.kin.thetaRangesDeg); ...
        'kin.kneeSlot',  isequal(A.kin.kneeSlot, B.kin.kneeSlot); ...
        'kin.pointFolds', isequal(A.kin.pointFolds, B.kin.pointFolds); ...
        'kin.kneeFit',   isequalt(A.kin.kneeFit, B.kin.kneeFit)};
    for c = 1:size(checks, 1)
        label = checks{c, 1};
        val = checks{c, 2};
        if ~val
            nDiff = nDiff + 1;
            fprintf('DIFF %s : %s\n', A.name, label);
        end
    end
    fprintf('%-12s %s\n', A.name, iIf(nDiff == 0, 'identical (all fields)', ...
        'see DIFF lines above'));
end
if nDiff > 0
    error('check_verbatim_20261023:mismatch', ...
        '%d field difference(s) across the 4 shared routes', nDiff)
end
fprintf('VERBATIM CHECK PASS: the 4 original routes are byte-identical.\n');
end

function tf = isequalt(a, b)
% Tolerant equality for the kneeFit struct (plain isequal also works; kept
% explicit so a future struct change cannot silently pass).
tf = isequal(a, b);
end

function s = iIf(cond, x, y)
if cond
    s = x;
else
    s = y;
end
end
