%EXPORT_SLX_TO_R2025A  Export SNS_Simscape models as R2025a-format copies.
% Creates <name>_R2025a.slx next to the originals (originals untouched).
% KneeReflexDemo's library links are re-pointed to the exported library copy,
% because an exported slx is internally renamed to match its file name.
% Run:  cd to this folder, then  matlab -batch "export_slx_to_R2025a"
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);
addpath(thisDir);

% 1) Library copy first (the demo relink depends on it).
load_system('SNS_Library.slx');
if strcmpi(get_param('SNS_Library', 'Lock'), 'on')
    set_param('SNS_Library', 'Lock', 'off');
end
Simulink.exportToVersion('SNS_Library', 'SNS_Library_R2025a', 'R2025a');
close_system('SNS_Library', 0);
fprintf('EXPORTED_OK: SNS_Library_R2025a.slx\n');

% 2) Demo: relink to the renamed library copy, then export.
load_system('SNS_Library_R2025a');   % keep loaded so the relink resolves
load_system('KneeReflexDemo.slx');
blks = find_system('KneeReflexDemo', 'LookUnderMasks', 'all', 'Type', 'Block');
nLink = 0; nRelink = 0; nFail = 0;
for h = 1:numel(blks)
    b = blks{h};
    try
        if ~strcmp(get_param(b, 'LinkStatus'), 'none')
            nLink = nLink + 1;
            rb = get_param(b, 'ReferenceBlock');
            if strncmp(rb, 'SNS_Library/', numel('SNS_Library/'))
                set_param(b, 'ReferenceBlock', ['SNS_Library_R2025a/' rb(numel('SNS_Library/')+1:end)]);
                nRelink = nRelink + 1;
            end
        end
    catch ME
        nFail = nFail + 1;
        fprintf('RELINK_FAILED: %s : %s\n', b, ME.message);
    end
end
fprintf('RELINKED %d of %d library links (%d failed)\n', nRelink, nLink, nFail);
Simulink.exportToVersion('KneeReflexDemo', 'KneeReflexDemo_R2025a', 'R2025a');
close_system('KneeReflexDemo', 0);
fprintf('EXPORTED_OK: KneeReflexDemo_R2025a.slx\n');

% 3) Any other plain models: export as-is.
others = dir('*.slx');
skip = {'SNS_Library.slx', 'KneeReflexDemo.slx'};
for k = 1:numel(others)
    [~, b] = fileparts(others(k).name);
    if any(strcmp(others(k).name, skip)) || endsWith(b, '_R2025a')
        continue;
    end
    try
        load_system(others(k).name);
        Simulink.exportToVersion(b, [b '_R2025a'], 'R2025a');
        fprintf('EXPORTED_OK: %s_R2025a.slx\n', b);
        close_system(b, 0);
    catch ME
        fprintf('FAILED: %s : %s\n', others(k).name, ME.message);
        try close_system(b, 0); catch; end
    end
end

% 4) Verify: reload the copies fresh; every demo link must resolve to the copy.
load_system('SNS_Library_R2025a');
load_system('KneeReflexDemo_R2025a');
blks = find_system('KneeReflexDemo_R2025a', 'LookUnderMasks', 'all', 'Type', 'Block');
nOK = 0; nBad = 0; nNone = 0;
for h = 1:numel(blks)
    ls = get_param(blks{h}, 'LinkStatus');
    switch ls
        case 'none'
            nNone = nNone + 1;
        case 'unresolved'
            nBad = nBad + 1;
            fprintf('BROKEN_LINK: %s\n', blks{h});
        otherwise
            if strncmp(get_param(blks{h}, 'ReferenceBlock'), 'SNS_Library_R2025a/', numel('SNS_Library_R2025a/'))
                nOK = nOK + 1;
            else
                nBad = nBad + 1;
                fprintf('WRONG_TARGET: %s -> %s\n', blks{h}, get_param(blks{h}, 'ReferenceBlock'));
            end
    end
end
fprintf('VERIFY: %d links OK, %d broken/wrong, %d unlinked blocks\n', nOK, nBad, nNone);
close_system('KneeReflexDemo_R2025a', 0);
close_system('SNS_Library_R2025a', 0);
fprintf('ALL_DONE\n');
