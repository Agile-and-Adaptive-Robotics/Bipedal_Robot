function probe_ankle_20260926()
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
sns = fileparts(fileparts(mfilename('fullpath')));
load_system(fullfile(sns, [mdl '.slx']));
for nm = {'ankle_R', 'foot_T_R', 'ankle_off_R', 'foot_R'}
    blk = [sub '/' nm{1}];
    fprintf('--- %s (exists=%d) ---\n', nm{1}, getSimulinkBlockHandle(blk) > 0);
    ph = get_param(blk, 'PortHandles');
    fn = fieldnames(ph);
    for k = 1:numel(fn)
        for j = 1:numel(ph.(fn{k}))
            h = ph.(fn{k})(j);
            ln = -1;
            try, ln = get_param(h, 'Line'); catch, end
            far = '<unconnected>';
            if ln > 0
                phs = [];
                try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
                br = -1;
                try, br = get_param(h, 'Branch'); catch, end
                if isscalar(br) && br > 0
                    try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
                end
                far = '';
                for q = 1:numel(phs)
                    hh = phs(q);
                    if isscalar(hh) && hh > 0 && hh ~= h
                        pn = '<rot?>';
                        try, pn = get_param(get_param(hh, 'Parent'), 'Name'); catch, end
                        far = [far ' ' pn]; %#ok<AGROW>
                    end
                end
                if isempty(far), far = '<branch-only>'; end
            end
            fprintf('  %s(%d) line=%d -> %s\n', fn{k}, j, ln, far);
        end
    end
end
close_system(mdl, 0);
end
