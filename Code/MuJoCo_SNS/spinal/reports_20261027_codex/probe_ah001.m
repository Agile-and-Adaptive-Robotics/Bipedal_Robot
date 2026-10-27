outDir = fileparts(mfilename('fullpath'));
repoRoot = outDir;
while ~strcmpi(string(lastFolder(repoRoot)), "Bipedal_Robot")
    parentDir = fileparts(repoRoot);
    assert(~strcmp(parentDir, repoRoot), 'Repository root not found.');
    repoRoot = parentDir;
end
diary(fullfile(outDir, 'ah001_matlab_probe.log'));
disp(version);
modelDir = fullfile(repoRoot, 'Code', 'Matlab', 'SNS_Simscape');
addpath(modelDir);
source = fullfile(modelDir, 'mdl_humanoid_lower_ah001_imported_R2025a.slx');
load_system(source);
mdl = 'mdl_humanoid_lower_ah001_imported_R2025a';
report = struct;
report.source = source;
report.release = version('-release');
report.initFcn = get_param(mdl, 'InitFcn');
report.preLoadFcn = get_param(mdl, 'PreLoadFcn');
report.sourceBlocks = find_system(mdl, 'LookUnderMasks', 'all', ...
    'FollowLinks', 'off', 'BlockType', 'Reference');
report.sourceReferences = cellfun(@(b) get_param(b, 'SourceBlock'), ...
    report.sourceBlocks, 'UniformOutput', false);
try
    load_system('sm_lib');
    set_param(mdl, 'SimulationCommand', 'update');
    report.compilePassed = true;
    report.compileError = '';
catch err
    report.compilePassed = false;
    report.compileError = getReport(err, 'extended', 'hyperlinks', 'off');
    disp(report.compileError);
end
fid = fopen(fullfile(outDir, 'ah001_matlab_probe.json'), 'w');
fwrite(fid, jsonencode(report, 'PrettyPrint', true));
fclose(fid);
close_system(mdl, 0);
diary off;

function name = lastFolder(path)
[~, name] = fileparts(path);
end
