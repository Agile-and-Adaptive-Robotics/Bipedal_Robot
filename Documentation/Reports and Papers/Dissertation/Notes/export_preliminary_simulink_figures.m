% Export genuine saved Simulink diagrams; never simulate or save source models.
repoRoot = fileparts(mfilename('fullpath'));
while ~strcmp(getLastFolder(repoRoot), 'Bipedal_Robot')
    parent = fileparts(repoRoot);
    assert(~strcmp(parent, repoRoot), 'Repository root not found.');
    repoRoot = parent;
end
sourceDir = fullfile(repoRoot, 'Code', 'Matlab', 'SNS_Simscape');
outDir = fullfile(repoRoot, 'Documentation', 'Reports and Papers', ...
    'Dissertation', 'ProofFinal', 'figs', 'Preliminary');
addpath(sourceDir);
fprintf('MATLAB %s\n', version);
load_system(fullfile(sourceDir, 'SNS_Library.slx'));
load_system(fullfile(sourceDir, 'KneeReflexDemo.slx'));
cleanup = onCleanup(@closeWithoutSaving);
systems = {'SNS_Library', 'KneeReflexDemo', 'SNS_Library/NonSpikingNeuron'};
stems = {'simulink_sns_library', 'simulink_knee_reflex', 'simulink_neuron_detail'};
for k = 1:numel(systems)
    fprintf('Exporting %s\n', systems{k});
    if k == 3
        % Render native primitive icons inside the masked neuron before print.
        open_system(systems{k}, 'force');
        drawnow;
    end
    print(['-s' systems{k}], '-dpng', '-r300', fullfile(outDir, [stems{k} '.png']));
    print(['-s' systems{k}], '-dpdf', fullfile(outDir, [stems{k} '.pdf']));
end
blocks = find_system('SNS_Library', 'SearchDepth', 1, 'Type', 'Block');
fprintf('Saved library block count: %d\n', numel(blocks));
fprintf('%s\n', blocks{:});
% Native Simulink layout only: preserve all blocks, parameters, and connections.
beforeConnectivity = connectionSnapshot('KneeReflexDemo');
open_system('KneeReflexDemo');
Simulink.BlockDiagram.arrangeSystem('KneeReflexDemo', 'FullLayout', 'true');
drawnow;
assert(isequal(beforeConnectivity, connectionSnapshot('KneeReflexDemo')), ...
    'Diagram layout unexpectedly changed connectivity.');
print('-sKneeReflexDemo', '-dpng', '-r300', fullfile(outDir, 'simulink_knee_reflex_arranged.png'));
print('-sKneeReflexDemo', '-dpdf', fullfile(outDir, 'simulink_knee_reflex_arranged.pdf'));
fprintf('Reflex diagram auto-arranged in memory; connectivity verified unchanged.\n');
% The saved library has a very tall single column. Reposition the same seven
% masked blocks only in memory for a compact page figure, then discard changes.
set_param('SNS_Library', 'Lock', 'off');
orderedNames = {'NonSpikingNeuron', 'NonSpikingSynapse', 'SpikingLIFNeuron', ...
    'IaMuscleSpindle', 'IbGolgiTendon', 'MuscleActivation', 'BPAForce'};
for k = 1:numel(orderedNames)
    col = mod(k-1, 3); row = floor((k-1)/3);
    x = 40 + 220*col; y = 40 + 145*row;
    set_param(['SNS_Library/' orderedNames{k}], 'Position', [x y x+160 y+90]);
end
set_param('SNS_Library', 'Lock', 'on');
print('-sSNS_Library', '-dpng', '-r300', fullfile(outDir, 'simulink_sns_library_grid.png'));
print('-sSNS_Library', '-dpdf', fullfile(outDir, 'simulink_sns_library_grid.pdf'));
fprintf('No simulation or model saves requested. Grid and reflex arrangement change display layout only, in memory.\n');
clear cleanup;

function leaf = getLastFolder(path)
    [~, leaf] = fileparts(path);
end

function closeWithoutSaving()
    close_system('KneeReflexDemo', 0);
    close_system('SNS_Library', 0);
end

function connections = connectionSnapshot(model)
    names = sort(find_system(model, 'SearchDepth', 1, 'Type', 'Block'));
    connections = {};
    for n = 1:numel(names)
        handles = get_param(names{n}, 'PortHandles');
        for p = 1:numel(handles.Inport)
            line = get_param(handles.Inport(p), 'Line');
            sourceName = '<unconnected>'; sourcePort = -1;
            if line ~= -1
                source = get_param(line, 'SrcPortHandle');
                if source ~= -1
                    sourceName = get_param(source, 'Parent');
                    sourcePort = get_param(source, 'PortNumber');
                end
            end
            connections(end+1,:) = {names{n}, p, sourceName, sourcePort}; %#ok<AGROW>
        end
    end
end
