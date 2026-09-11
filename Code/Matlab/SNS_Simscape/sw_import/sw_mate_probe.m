% sw_mate_probe.m — MATLAB COM probe of 09_BA_003 mates (ActiveDoc must be
% the assembly; MATLAB handles SW byref GetEntity args natively). Read-only.
asmTitle = '09_BA_003.SLDASM';
sw = actxserver('SLDWorks.Application');
try, sw.Visible = 0; catch, end
doc = invoke(sw, 'ActiveDoc');
if isempty(doc)
    error('Open the assembly in SolidWorks first');
end
fprintf('active: %s\n', invoke(doc, 'GetTitle'));

comps = invoke(doc, 'GetComponents', false);
fprintf('components: %d\n', numel(comps));

seen = {};
nMates = 0;
for ci = 1:numel(comps)
    c = comps{ci};     % GetComponents returns a cell array
    cname = invoke(c, 'Name2');
    ml = invoke(c, 'GetMates');
    if isempty(ml), continue; end
    for mi = 1:numel(ml)
        m = ml{mi};
        mtype = invoke(m, 'Type');
        nMates = nMates + 1;
        for side = 0:1
            ocName = '?'; et = -1; ok = false;
            try
                e = invoke(m, 'GetEntity', int32(side));
                if ~isempty(e)
                    oc = e.GetComponent;
                    if ~isempty(oc), ocName = oc.Name2; end
                    et = e.GetType;
                    ok = true;
                end
            catch ME
                ocName = ['ERR: ' ME.message(1:min(60, end))];
            end
            key = sprintf('%d|%s|%s|%d|%d', mtype, cname, ocName, et, side);
            if ~any(strcmp(seen, key))
                seen{end+1} = key; %#ok<*AGROW>
                fprintf('mate[type %d] on %-24s side %d: comp=%s selType=%d\n', ...
                    mtype, cname, side, ocName, et);
            end
        end
    end
end
fprintf('mate instances seen: %d, unique rows: %d\n', nMates, numel(seen));
fprintf('MATE PROBE DONE\n');
