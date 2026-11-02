% mass_ledger.m  (AH001_Simscape_Repair, 2026-11-01)
% Sum the 75 File Solid instance masses in the repaired model, group by part
% ID, and compare against the old record (5.3566638205 kg incl. feet).

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'mass_ledger_20261101.txt');
fid = fopen(logFile, 'w');

mdl = 'ah001_repair';
load_system(fullfile(repairRoot, 'working', [mdl '.slx']));

hws = get_param(mdl, 'ModelWorkspace');
sm = evalin(hws, 'smiData');

blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try, ref(k) = strtrim(strrep(string(get_param(blks{k}, 'ReferenceBlock')), newline, " ")); catch, ref(k) = ""; end
end
fs = blks(ref == "sm_lib/Body Elements/File Solid");

p(fid, 'File Solid instances: %d\n', numel(fs));
total = 0;
byPart = containers.Map('KeyType','char','ValueType','double');
rows = {};
for k = 1:numel(fs)
    % Mass param is an expression like smiData.Solid(4).mass -> resolve via ws
    mexpr = get_param(fs{k}, 'Mass');
    mi = regexp(mexpr, 'Solid\((\d+)\)', 'tokens', 'once');
    if isempty(mi)
        m = str2double(mexpr);
        id = 'literal';
    else
        i = str2double(mi{1});
        m = sm.Solid(i).mass;
        try, id = char(string(sm.Solid(i).ID)); catch, id = sprintf('Solid(%d)', i); end
    end
    total = total + m;
    if isKey(byPart, id)
        byPart(id) = byPart(id) + m;
    else
        byPart(id) = m;
    end
    rows{end+1} = sprintf('%12.9f kg  %-64s %s', m, id, strrep(blks{k}, [mdl '/'], '')); %#ok<AGROW>
end
% note: instance count per part
cnt = containers.Map('KeyType','char','ValueType','double');
for k = 1:numel(rows)
    tok = regexp(rows{k}, 'kg\s+(\S+)', 'tokens', 'once');
    if ~isempty(tok)
        key = tok{1};
        if isKey(cnt, key), cnt(key) = cnt(key) + 1; else, cnt(key) = 1; end
    end
end

p(fid, '\nper-instance ledger:\n');
for k = 1:numel(rows), p(fid, '  %s\n', rows{k}); end
p(fid, '\nper-part totals:\n');
ids = byPart.keys;
for k = 1:numel(ids)
    key = ids{k};
    n = 0;
    for j = 1:numel(rows)
        if contains(rows{j}, key), n = n + 1; end
    end
    p(fid, '  %12.9f kg  x%d  %s\n', byPart(key), n, key);
end
p(fid, '\nTOTAL lower-body mass (CAD instances, NO feet, NO spine/payload): %.10f kg\n', total);
p(fid, 'old record incl. feet: 5.3566638205 kg (different assembly revision - material changes)\n');
fclose(fid);
close_system(mdl, 0);

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
