% Timing: MuJoCo stepped from Simulink (bridge S-function, the shipped
% mj_gettingStarted example) vs pure-Python MuJoCo 2.3.7, same 2 s @ 1 ms
% on the same small model (the example's cartpole MJCF if present, else
% the humanoid shipped with the blockset's mujoco lib).
here = fileparts(mfilename('fullpath'));
root = 'D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\mujoco_bridge\blockset\mujoco-simulink-blockset-main';
cd(root); addpath(root); addpath(fullfile(root, 'blocks'));

mdl = 'mj_gettingStarted';
load_system(fullfile(root, 'examples', [mdl '.slx']));
xml = get_param(find_system(mdl, 'Mask', 'on'), 'MaskValues'); %#ok<NASGU>
% find the mj block and read its xml + sample time
mjb = find_system(mdl, 'LookUnderMasks', 'all', 'ReferenceBlock', '');
t0 = tic;
out = sim(mdl, 'StopTime', '2', 'ReturnWorkspaceOutputs', 'on');
tSim = toc(t0);
close_system(mdl, 0);
fprintf('SIMULINK-BRIDGE: 2 s sim wall time %.3f s\n', tSim);
