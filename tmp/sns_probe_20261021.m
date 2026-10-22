% Probe: census of SNS_Library.slx as saved (blocks, masks, annotations, layout).
here = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape';
addpath(here);
fprintf('MATLAB %s\n', version);
load_system(fullfile(here, 'SNS_Library.slx'));
blocks = find_system('SNS_Library', 'SearchDepth', 1, 'Type', 'Block');
fprintf('BLOCK COUNT: %d\n', numel(blocks));
for k = 1:numel(blocks)
    nm = get_param(blocks{k}, 'Name');
    try, mt = get_param(blocks{k}, 'MaskType'); catch, mt = ''; end
    p = get_param(blocks{k}, 'Position');
    try, mps = strcat(get_param(blocks{k}, 'MaskNames')', '=');
    catch, mps = {}; end
    fprintf('  %-22s | mask=%-28s | pos=[%d %d %d %d]\n', nm, mt, p(1), p(2), p(3), p(4));
    if ~isempty(mps), fprintf('      params: %s\n', strjoin(mps, ' ')); end
end
anns = find_system('SNS_Library', 'SearchDepth', 1, 'Type', 'annotation');
fprintf('ANNOTATION COUNT: %d\n', numel(anns));
for k = 1:numel(anns)
    t = get_param(anns{k}, 'Text');
    p = get_param(anns{k}, 'Position');
    fprintf('  ann pos=[%d %d] : %s\n', p(1), p(2), strrep(t, newline, ' | '));
end
fprintf('LOCKED: %s\n', get_param('SNS_Library', 'Lock'));
close_system('SNS_Library', 0);
fprintf('PROBE DONE\n');
