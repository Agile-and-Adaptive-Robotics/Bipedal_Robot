% smimport_foot.m  (AH001_Simscape_Repair, 2026-11-01)
% Import the Foot Ankle CAD-XML and census what the ankle/MTP mates became.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
expDir = fullfile(repairRoot, 'cad', 'export_foot');
logFile = fullfile(repairRoot, 'logs', 'smimport_foot_20261101.txt');
fid = fopen(logFile, 'w');

src = fullfile(expDir, 'foot_ankle.xml.xml');
dst = fullfile(expDir, 'foot_ankle.xml');
if exist(src, 'file') && ~exist(dst, 'file'), copyfile(src, dst); end

cd(expDir);
mdlName = 'foot_import';
smimport('foot_ankle.xml', 'ModelName', mdlName);
p(fid, 'imported\n');

blks = find_system(mdlName, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try, ref(k) = strtrim(strrep(string(get_param(blks{k}, 'ReferenceBlock')), newline, " ")); catch, ref(k) = ""; end
end
p(fid, 'total blocks: %d\n', numel(blks));
uj = unique(ref(startsWith(ref, "sm_lib/")));
for t = uj'
    p(fid, '  %2d x %s\n', nnz(ref == t), t);
end

% masses from the generated data file
df = dir(fullfile(expDir, '*_DataFile.m'));
p(fid, 'data files: %s\n', strjoin({df.name}, ', '));
try
    sm = evalin('caller', '');  %#ok<NASGU>
catch
end

save_system(mdlName, fullfile(repairRoot, 'working', [mdlName '.slx']));
close_system(mdlName, 0);
p(fid, 'FOOT IMPORT DONE\n');
fclose(fid);

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
