% probe2_compile_topo.m  (AH001_Simscape_Repair, 2026-11-01)
% With SNS_Library on the path: resolve link status, compare saved model-workspace
% smiData vs DataFile smiData, dump joint/force topology, attempt compile with a
% full extended error report. Also census the R2025a export variant.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
workingDir = fullfile(repairRoot, 'working');
snsDir     = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape';
logFile    = fullfile(repairRoot, 'logs', 'probe2_compile_topo_20261101.txt');
fid = fopen(logFile, 'w');

addpath(snsDir);

%% Main model
mdl = 'mdl_humanoid_lower_ah001_imported';
load_system(fullfile(workingDir, [mdl '.slx']));
p(fid, 'loaded %s\n', mdl);

% --- saved smiData vs data file smiData
hws = get_param(mdl, 'ModelWorkspace');
savedOK = false;
try
    smSaved = evalin(hws, 'smiData');
    savedOK = true;
    p(fid, 'saved model-workspace smiData present\n');
catch
    p(fid, 'saved model-workspace smiData ABSENT\n');
end
evalin(hws, sprintf('clear; run(''%s'');', fullfile(workingDir, 'mdl_humanoid_lower_ah001_DataFile.m')));
smFile = evalin(hws, 'smiData');
if savedOK
    same = isequal(smSaved, smFile);
    p(fid, 'saved smiData == DataFile smiData : %d\n', same);
end

% --- link status of every block
blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
p(fid, '\n--- unresolved links ---\n');
nUnres = 0;
for k = 1:numel(blks)
    try
        ls = get_param(blks{k}, 'LinkStatus');
        if ~any(strcmp(ls, {'none','implicit'}))
            nUnres = nUnres + 1;
            p(fid, '  %s  [%s]\n', blks{k}, ls);
        end
    catch
    end
end
p(fid, 'unresolved-link blocks: %d\n', nUnres);

% --- joint topology dump
p(fid, '\n--- joint blocks ---\n');
jointRefPat = "sm_lib/Joints/";
for k = 1:numel(blks)
    try
        r = string(get_param(blks{k}, 'ReferenceBlock'));
    catch
        continue;
    end
    if startsWith(r, jointRefPat) || contains(r, "Joint")
        nm = get_param(blks{k}, 'Name');
        parent = fileparts(blks{k});
        try
            ndof = get_param(blks{k}, 'NumberDoF');
        catch
            ndof = '?';
        end
        p(fid, '  %-28s in %-42s DoF=%s  ref=%s\n', nm, parent, ndof, r);
    end
end

% --- internal force + sensor blocks
p(fid, '\n--- force/sensor blocks ---\n');
for k = 1:numel(blks)
    try
        r = string(get_param(blks{k}, 'ReferenceBlock'));
    catch
        continue;
    end
    if contains(r, "Internal Force") || contains(r, "Transform Sensor")
        p(fid, '  %s\n', blks{k});
    end
end

% --- solids: names + masses from smiData.Solid
p(fid, '\n--- smiData.Solid census (mass in source units) ---\n');
for k = 1:numel(smFile.Solid)
    s = smFile.Solid(k);
    nm = '';
    try, nm = char(s.ID); catch, end
    m = NaN; try, m = s.mass; catch, end
    p(fid, '  %2d mass=%-12.6g ID=%s\n', k, m, nm);
end

%% Compile attempt with full report
p(fid, '\n=== sim(StopTime=0) attempt (links resolved) ===\n');
try
    simOut = sim(mdl, 'StopTime', '0'); %#ok<NASGU>
    p(fid, 'sim: OK\n');
catch ME
    rep = getReport(ME, 'extended', 'hyperlinks', 'off');
    p(fid, 'sim FAIL:\n%s\n', rep);
end
try, close_system(mdl, 0); catch, end

%% R2025a variant census
p(fid, '\n=== R2025a variant (mdl_humanoid_lower_ah001_imported_R2025a) ===\n');
mdl2 = 'mdl_humanoid_lower_ah001_imported_R2025a';
try
    load_system(fullfile(workingDir, [mdl2 '.slx']));
    cbn = {'PreLoadFcn','InitFcn','PostLoadFcn'};
    for k = 1:numel(cbn)
        p(fid, '%-12s = [%s]\n', cbn{k}, strtrim(get_param(mdl2, cbn{k})));
    end
    hws2 = get_param(mdl2, 'ModelWorkspace');
    v2 = whos(hws2);
    p(fid, 'R2025a model workspace vars: %d\n', numel(v2));
    b2 = find_system(mdl2, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
    p(fid, 'R2025a total blocks: %d\n', numel(b2));
    r2 = strings(numel(b2),1);
    for k = 1:numel(b2)
        try, r2(k) = string(get_param(b2{k}, 'ReferenceBlock')); catch, r2(k) = ""; end
    end
    [u2,~,i2] = unique(r2);
    for k = 1:numel(u2)
        if startsWith(u2(k), "sm_lib/") || contains(u2(k), "SNS")
            p(fid, '  %3d x %s\n', nnz(i2==k), u2(k));
        end
    end
    try, close_system(mdl2, 0); catch, end
catch ME
    p(fid, 'R2025a census FAIL: %s\n', ME.message);
end

p(fid, '\nPROBE2 DONE.\n');
fclose(fid);

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
