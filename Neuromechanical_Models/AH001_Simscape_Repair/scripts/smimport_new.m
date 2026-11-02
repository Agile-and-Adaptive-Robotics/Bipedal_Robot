% smimport_new.m  (AH001_Simscape_Repair, 2026-11-01)
% Import the fresh CAD-XML export (from the real 10_AH_001.SLDASM) into a new
% Simscape Multibody model, save it into working/, and census the result.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
expDir = fullfile(repairRoot, 'cad', 'export');
logFile = fullfile(repairRoot, 'logs', 'smimport_new_20261101.txt');
fid = fopen(logFile, 'w');

% normalize the doubled extension the add-in produced
src = fullfile(expDir, 'ah001_full.xml.xml');
dst = fullfile(expDir, 'ah001_full.xml');
if exist(src, 'file') && ~exist(dst, 'file')
    copyfile(src, dst);
end
p(fid, 'xml ready: %d\n', exist(dst, 'file') == 2);

cd(expDir);   % STEP paths are relative; import must run in the export folder

mdlName = 'ah001_full_import';
p(fid, 'running smimport ...\n');
try
    smimport('ah001_full.xml', 'ModelName', mdlName);
    p(fid, 'smimport returned without error\n');
catch ME
    p(fid, 'smimport ERR: %s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    fclose(fid);
    rethrow(ME);
end

% find the created model (name per skill notes can differ)
cands = find_system('type', 'block_diagram');
p(fid, 'open models: %s\n', strjoin(cands, ', '));
mdl = mdlName;
if ~any(strcmp(cands, mdl))
    p(fid, 'WARNING: expected model %s not in list\n', mdl);
end

% census
blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
p(fid, 'total blocks: %d\n', numel(blks));
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try, ref(k) = strtrim(string(get_param(blks{k}, 'ReferenceBlock'))); catch, ref(k) = ""; end
end
[uref,~,iu] = unique(ref);
for k = 1:numel(uref)
    if startsWith(uref(k), "sm_lib/")
        p(fid, '  %3d x %s\n', nnz(iu==k), uref(k));
    end
end

% save into working/
outPath = fullfile(repairRoot, 'working', [mdlName '.slx']);
try
    save_system(mdl, outPath);
    p(fid, 'saved %s\n', outPath);
catch ME
    p(fid, 'save ERR: %s\n', ME.message);
end

try, close_system(mdl, 0); catch, end
p(fid, 'IMPORT DONE.\n');
fclose(fid);

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
