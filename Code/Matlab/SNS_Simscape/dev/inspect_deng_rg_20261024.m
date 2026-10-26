% List the actual blocks of the freshly built SNS_Deng_RG (read-only).
root = 'D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape';
cd(root); addpath(root);
mdl = 'SNS_Deng_RG';
load_system(fullfile(root, 'demos', [mdl '.slx']));
blks = find_system(mdl, 'SearchDepth', 1, 'Type', 'Block');
fprintf('--- %s top-level blocks (%d) ---\n', mdl, numel(blks));
for k = 1:numel(blks)
    fprintf('%s | %s\n', blks{k}, get_param(blks{k}, 'BlockType'));
end
% HC port count check
hc = [mdl '/pair1_HC_ext'];
ph = get_param(hc, 'PortHandles');
fprintf('HC_ext: %d inports, %d outports\n', numel(ph.Inport), numel(ph.Outport));
close_system(mdl, 0);
disp('inspect done');
