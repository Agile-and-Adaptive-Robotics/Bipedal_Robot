% probe_frames_20261027.m -- frame reconciliation probe (read-only)
R = 'D:\GitHub\Bipedal_Robot';
D = fullfile(R, 'tmp', 'old_code_copy_20261022');
M = load(fullfile(D, 'MovePoints.mat'), 'knee_angle_x', 'knee_x', 'knee_angle_y', ...
    'knee_y', 'rect_fem_x', 'rect_fem_y', 'rect_fem_xD', 'rect_fem_yD');
fprintf('MovePoints knee_angle_x (rad):\n'); disp(M.knee_angle_x(:).');
fprintf('MovePoints knee_x (m):\n'); disp(M.knee_x(:).');
fprintf('MovePoints knee_angle_y (rad):\n'); disp(M.knee_angle_y(:).');
fprintf('MovePoints knee_y (m):\n'); disp(M.knee_y(:).');
S = load(fullfile(D, 'AttachPoints_Robot.mat'));
for nm = {'t1','t2','t3','t4','t5','t6','f1','f2','twr1','twr2','twr3','p10'}
    fprintf('%-6s = [%s]\n', nm{1}, num2str(S.(nm{1}).', 6));
end
fprintf('Patella-related vars in AttachPoints: ');
fn = fieldnames(S); disp(fn(regexpi(fn, 'pa')));
fprintf('FRAMES PROBE DONE\n');
