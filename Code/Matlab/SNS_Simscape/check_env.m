% Environment check: toolboxes + smlink presence
fprintf('MATLAB %s\n', version);
prods = {ver().Name};
want = {'Simulink','Simscape','Simscape Multibody','Simscape Multibody Link','Simscape Multibody Link for SolidWorks','Control System Toolbox','MATLAB Coder'};
for i = 1:numel(want)
    fprintf('%-40s installed=%d license=%d\n', want{i}, any(strcmpi(prods, want{i})), license('test', strrep(want{i},' ','_')));
end
smlinkdir = fullfile(matlabroot,'toolbox','physmod','smlink');
fprintf('smlink dir exists: %d  (%s)\n', isfolder(smlinkdir), smlinkdir);
if isfolder(fullfile(smlinkdir,'mw'))
    d = dir(fullfile(smlinkdir,'mw','*.m'));
    for k=1:numel(d), fprintf('  mw/%s\n', d(k).name); end
end
