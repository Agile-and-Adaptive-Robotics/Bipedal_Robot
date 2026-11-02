% probe_env_model.m  (AH001_Simscape_Repair, laptop DESKTOP-5Q16KE9, 2026-11-01)
% Gate-0/1 probe: MATLAB + Simscape Multibody runtime license, then inspect the
% copied import model's callbacks / model workspace / block census, run the
% DataFile in the model-workspace context, and attempt one compile.
% Read-only with respect to the ORIGINAL repo sources; operates on the copies
% in AH001_Simscape_Repair/working.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
workingDir = fullfile(repairRoot, 'working');
logFile    = fullfile(repairRoot, 'logs', 'probe_env_model_20261101.txt');
fid = fopen(logFile, 'w');

%% 1) Environment
p(fid, '=== MATLAB environment ===\n');
p(fid, 'MATLAB version: %s\n', version);
v = ver;
for k = 1:numel(v)
    nm = string(v(k).Name);
    if contains(nm, "Simscape", 'IgnoreCase', true) || contains(nm, "SimMechanics", 'IgnoreCase', true)
        p(fid, '  product: %s  version %s\n', v(k).Name, v(k).Version);
    end
end

%% 2) Runtime license probes (license('test') is unreliable per skill/handoff)
p(fid, '\n=== Simscape Multibody runtime probes ===\n');
try
    load_system('sm_lib');
    p(fid, 'load_system(''sm_lib''): OK\n');
catch ME
    p(fid, 'load_system(''sm_lib'') FAIL: %s\n', ME.message);
end
scratch = 'ah001_scratch_probe';
try
    new_system(scratch);
    add_block('sm_lib/Utilities/Mechanism Configuration', [scratch '/MechCfg']);
    p(fid, 'scratch Mechanism Configuration add_block: OK (license usable at runtime)\n');
catch ME
    p(fid, 'scratch add_block FAIL: %s\n', ME.message);
end
try
    close_system(scratch, 0);
catch
end

%% 3) Load the copied model and inspect callbacks + model workspace
p(fid, '\n=== Model callbacks / workspace (copy of mdl_humanoid_lower_ah001_imported.slx) ===\n');
mdl = 'mdl_humanoid_lower_ah001_imported';
try
    load_system(fullfile(workingDir, [mdl '.slx']));
catch ME
    p(fid, 'load_system FAIL: %s\n', ME.message);
    fclose(fid);
    rethrow(ME);
end
cbNames = {'PreLoadFcn','PostLoadFcn','InitFcn','StartFcn','StopFcn'};
for k = 1:numel(cbNames)
    cb = strtrim(get_param(mdl, cbNames{k}));
    p(fid, '%-12s = [%s]\n', cbNames{k}, cb);
end

hws = get_param(mdl, 'ModelWorkspace');
vars = whos(hws);
p(fid, 'model workspace variables before DataFile: %d\n', numel(vars));
for k = 1:numel(vars)
    p(fid, '  %-10s %s\n', vars(k).class, vars(k).name);
end

p(fid, 'model Dirty: %s\n', get_param(mdl, 'Dirty'));

%% 4) Block census (all reference-block histogram)
p(fid, '\n=== Block census ===\n');
blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
p(fid, 'total blocks: %d\n', numel(blks));
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try
        ref(k) = string(get_param(blks{k}, 'ReferenceBlock'));
    catch
        ref(k) = "";  % non-library block
    end
end
[uref, ~, iu] = unique(ref);
p(fid, 'reference-block histogram (%d distinct):\n', numel(uref));
for k = 1:numel(uref)
    p(fid, '  %3d x %s\n', nnz(iu == k), uref(k));
end

%% 5) Run the DataFile in the model workspace, then re-check
p(fid, '\n=== DataFile execution in model workspace ===\n');
dataFile = fullfile(workingDir, 'mdl_humanoid_lower_ah001_DataFile.m');
try
    evalin(hws, sprintf('clear; run(''%s'');', dataFile));
    vars2 = whos(hws);
    p(fid, 'DataFile run: OK; workspace variables after: %d\n', numel(vars2));
    for k = 1:numel(vars2)
        p(fid, '  %-10s %s\n', vars2(k).class, vars2(k).name);
    end
    sm = evalin(hws, 'smiData');
    fn = fieldnames(sm);
    p(fid, 'smiData fields: %s\n', strjoin(fn', ', '));
    for k = 1:numel(fn)
        try
            n = numel(sm.(fn{k}));
            p(fid, '  smiData.%s : %d entries\n', fn{k}, n);
        catch
            p(fid, '  smiData.%s : (scalar)\n', fn{k});
        end
    end
catch ME
    p(fid, 'DataFile run FAIL: %s\n', ME.message);
end

%% 6) One compile attempt (sim with zero stop time) on the copy
p(fid, '\n=== Compile/sim(StopTime=0) attempt ===\n');
try
    simOut = sim(mdl, 'StopTime', '0'); %#ok<NASGU>
    p(fid, 'sim(StopTime=0): OK\n');
catch ME
    p(fid, 'sim FAIL: %s\n', ME.message);
    for k = 1:min(numel(ME.stack), 6)
        p(fid, '  at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
    end
    try
        p(fid, 'cause: %s\n', ME.cause.message);
    catch
    end
end

try
    close_system(mdl, 0);
catch
end
p(fid, '\nPROBE DONE.\n');
fclose(fid);

%% Local function (must be last in a script file)
function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);        % stdout
    fprintf(fid, '%s', s);   % log file
end
