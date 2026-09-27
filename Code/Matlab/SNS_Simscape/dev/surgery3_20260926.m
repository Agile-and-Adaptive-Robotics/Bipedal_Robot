function surgery3_20260926()
%% surgery3_20260926  Seal the humanoid subsystem:
%%   1. root level: delete ALL loose hardware bodies, their joints, root
%%      Transforms, and the root World/Mechanism/Solver config blocks
%%      (keep only the x09_BA_001_1 subsystem)
%%   2. inside sub: delete leftover F#/Transform# frame blocks (world ties)
%%   3. add Solver Configuration + Mechanism Configuration (gravity -y) inside
%%      sub, branched onto the World_pelvis net
%%   4. compile + save

sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));
set_param(mdl, 'Lock', 'off');

%% ---- 1. clear root level ----------------------------------------------------
root = find_system(mdl, 'SearchDepth', 1, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
keep = {sub};
nDel = 0;for k = 1:numel(root)
    if any(strcmp(root{k}, keep)), continue; end
    % guard: never delete the model's solver/mech if named differently (we
    % re-add them inside sub anyway, so deleting here is fine)
    try
        delete_block(root{k});
        nDel = nDel + 1;
    catch ME
        fprintf('root delete failed %s: %s\n', root{k}, ME.message);
    end
end
fprintf('root level: deleted %d blocks (kept only the biped subsystem)\n', nDel);

%% ---- 2. clear leftover frame blocks inside sub ------------------------------
inner = find_system(sub, 'SearchDepth', 1, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
nDel2 = 0;
for k = 1:numel(inner)
    bn = get_param(inner{k}, 'Name');
    if ~isempty(regexpi(bn, '^(F\d+|Transform\d*)$', 'once'))
        try
            delete_block(inner{k});
            nDel2 = nDel2 + 1;
        catch ME
            fprintf('inner delete failed %s: %s\n', bn, ME.message);
        end
    end
end
fprintf('inside sub: deleted %d F#/Transform frame blocks\n', nDel2);

%% ---- 3. solver + mechanism inside sub ---------------------------------------
if isempty(find_system(sub, 'SearchDepth', 1, 'MaskType', 'Solver Configuration'))
    add_block('nesl_utility/Solver Configuration', [sub '/Solver Configuration'], ...
        'Position', [40 480 100 540]);
end
if isempty(find_system(sub, 'SearchDepth', 1, 'MaskType', 'Mechanism Configuration'))
    add_block('sm_lib/Utilities/Mechanism Configuration', [sub '/Mechanism Configuration'], ...
        'Position', [40 560 120 610]);
    set_param([sub '/Mechanism Configuration'], 'GravityVector', '[0 -9.80665 0]');
end
wiph = get_param([sub '/World_pelvis'], 'PortHandles');
wiPort = [wiph.RConn wiph.LConn];
sol = get_param([sub '/Solver Configuration'], 'PortHandles');
solP = [sol.RConn sol.LConn];
mec = get_param([sub '/Mechanism Configuration'], 'PortHandles');
mecP = [mec.RConn mec.LConn];
add_line(sub, solP(1), wiPort(1));   % branch onto the world net
add_line(sub, mecP(1), wiPort(1));   % branch onto the world net
fprintf('solver + mechanism configuration added inside sub (gravity -y)\n');

%% ---- 4. compile + save ------------------------------------------------------
try
    set_param(mdl, 'StopTime', '0.01');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('%s: UPDATE OK\n', mdl);
catch ME
    fprintf('%s: UPDATE FAILED: %s\n', mdl, ME.message);
    for c = 1:min(numel(ME.cause), 10)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 180)));
    end
end
save_system(mdl);
close_system(mdl, 0);
fprintf('=== surgery3 DONE ===\n');
end
