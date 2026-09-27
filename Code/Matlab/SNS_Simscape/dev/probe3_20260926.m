function probe3_20260926()
% replicate surgery2 section 3 verbatim with verbose prints
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

bodies = find_system(sub, 'SearchDepth', 1, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
fprintf('bodies found: %d\n', numel(bodies));
for k = 1:numel(bodies)
    bn = get_param(bodies{k}, 'Name');
    if numel(bn) < 6 || ~strcmpi(bn(end-4:end), '_RIGID'), continue; end
    ph = get_param(bodies{k}, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    fprintf('BODY %-46s ports=%d\n', bn, numel(ports));
    if strcmp(bn, 'x02_01_PE_001_1_RIGID')
        for pp = 1:min(numel(ports), 8)
            far = far_ends(ports(pp));
            fprintf('  port%d far_ends rows=%d', pp, size(far, 1));
            for f = 1:size(far, 1)
                fprintf(' [%s]', far{f, 1});
            end
            fprintf('\n');
            if size(far, 1) == 0
                ln = get_param(ports(pp), 'Line');
                fprintf('    (line=%.0f -> far_ends found nothing)\n', ln);
            end
        end
    end
end
close_system(mdl, 0);
end

function far = far_ends(porth)
far = {};
ln = -1;
try, ln = get_param(porth, 'Line'); catch, end
if ln < 0, return; end
phs = [];
try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
br = -1;
try, br = get_param(porth, 'Branch'); catch, end
if isscalar(br) && br > 0
    try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
end
seen = [];
for q = 1:numel(phs)
    h = phs(q);
    if ~isscalar(h) || h <= 0 || h == porth, continue; end
    if ~isempty(seen) && any(seen == h), continue; end
    seen(end+1) = h; %#ok<AGROW>
    try
        bn = get_param(get_param(h, 'Parent'), 'Name');
    catch
        continue;
    end
    far(end+1, :) = {bn, h}; %#ok<AGROW>
end
end
