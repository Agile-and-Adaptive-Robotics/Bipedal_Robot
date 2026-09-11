% dump_wiring2.m — authoritative line-by-line dump of BioMuscle + BPA_10mm
load_system('SNS_Library');
for sys = {'SNS_Library/BioMuscle', 'SNS_Library/BPA_10mm'}
    fprintf('==== %s ====\n', sys{1});
    lines = find_system(sys{1}, 'SearchDepth', 1, 'LookUnderMasks', 'all', 'Type', 'line');
    for k = 1:numel(lines)
        lh = lines{k};
        try
            sp = get_param(lh, 'SrcPortHandle');
            if sp == -1, continue; end
            sb = get_param(sp, 'Parent');
            sn = regexprep(sb, [sys{1} '/'], '');
            spn = get_param(sp, 'PortNumber');
            dph = get_param(lh, 'DstPortHandle');
            for j = 1:numel(dph)
                db = get_param(dph(j), 'Parent');
                dn = regexprep(db, [sys{1} '/'], '');
                dpn = get_param(dph(j), 'PortNumber');
                fprintf('  %s/%d -> %s/%d\n', sn, spn, dn, dpn);
            end
        catch ME
            fprintf('  (line %d: %s)\n', k, ME.message);
        end
    end
end
% key params
for blk = {'SNS_Library/BioMuscle/concB', 'SNS_Library/BioMuscle/bAb', ...
           'SNS_Library/BioMuscle/hillA_c', 'SNS_Library/BioMuscle/lOpt_c', ...
           'SNS_Library/BioMuscle/Adeg', 'SNS_Library/BioMuscle/lslack', ...
           'SNS_Library/BioMuscle/vSwitch', 'SNS_Library/BioMuscle/concDiv'}
    fprintf('%s : %s\n', blk{1}, class(get_param(blk{1}, 'ObjectParameters')));
end
fprintf('concB Value=%s\n', get_param('SNS_Library/BioMuscle/concB', 'Value'));
fprintf('bAb   Value=%s\n', get_param('SNS_Library/BioMuscle/bAb', 'Value'));
fprintf('hillA Value=%s\n', get_param('SNS_Library/BioMuscle/hillA_c', 'Value'));
fprintf('vSw crit=%s\n', get_param('SNS_Library/BioMuscle/vSwitch', 'Criteria'));
