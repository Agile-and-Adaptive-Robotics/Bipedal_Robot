function probe_deng_dump()
% Dump the compiled contents of the linked HCNeuron inside SNS_Deng_RG:
% every block, its key parameter, and the connectivity.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
load_system('SNS_Deng_Library');
blk = 'SNS_Deng_Library/HCNeuron';
blks = find_system(blk, 'LookUnderMasks', 'all', 'Type', 'Block');
for i = 1:numel(blks)
    bn = strrep(blks{i}, [blk '/'], '');
    bt = get_param(blks{i}, 'BlockType');
    extra = '';
    try
        extra = [' Inputs=' get_param(blks{i}, 'Inputs')];
    catch
    end
    if isempty(extra)
        try
            extra = [' Value=' get_param(blks{i}, 'Value')];
        catch
        end
    end
    if isempty(extra)
        try
            extra = [' Gain=' get_param(blks{i}, 'Gain')];
        catch
        end
    end
    if isempty(extra)
        try
            extra = [' IC=' get_param(blks{i}, 'InitialCondition')];
        catch
        end
    end
    if isempty(extra)
        try
            extra = [' Op=' get_param(blks{i}, 'Operator')];
        catch
        end
    end
    fprintf('%-12s [%s]%s\n', bn, bt, extra);
end
fprintf('\nconnectivity:\n');
lh = find_system(blk, 'LookUnderMasks', 'all', 'Type', 'line');
for i = 1:numel(lh)
    sp = get_param(lh{i}, 'SrcPortHandle');
    dp = get_param(lh{i}, 'DstPortHandle');
    if ~isempty(sp) && ~isempty(dp)
        sb = get_param(sp, 'Parent');
        db = get_param(dp, 'Parent');
        fprintf('  %s -> %s\n', strrep(sb, [blk '/'], ''), ...
                strrep(db, [blk '/'], ''));
    end
end
close_system('SNS_Deng_RG', 0);
end
