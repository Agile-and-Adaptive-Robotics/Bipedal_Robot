% dump_wiring.m — print every block's input sources inside SNS_Library/BioMuscle
load_system('SNS_Library');
blks = find_system('SNS_Library/BioMuscle', 'SearchDepth', 1, 'LookUnderMasks', 'all', 'Type', 'Block');
for k = 1:numel(blks)
    nm = regexprep(blks{k}, 'SNS_Library/BioMuscle/', '');
    pc = get_param(blks{k}, 'PortConnectivity');
    for j = 1:numel(pc)
        p = pc(j);
        srcs = p.SrcBlock;
        if isempty(srcs) || srcs <= 0
            continue;
        end
        try
            sb = get_param(srcs, 'Parent');
            fprintf('%-12s <- %s\n', nm, regexprep(sb, 'SNS_Library/BioMuscle/', ''));
        catch
        end
    end
end
