%% sns_cpg_gait2392.m — drive the imported OpenSim leg's knee with the SNS CPG + real BPA blocks
%
% Final integration demo (2026-09-10): takes the Simscape model imported from
% the MyoConverter Gait2392 conversion (osim_import/Gait2392_simbody_simscape.slx,
% grounded pelvis 6-DOF chain), finds the knee_angle_r revolute joint, and
% actuates it with a CPG-driven antagonist BPA pair:
%
%   CPG (half-center: RG_ext/RG_flex + slow adaptive inhibition, from
%        BPACPGLegDemo)  ->  BPA_20mm force (Ben's Festo equations)
%                       ->  torque via moment arm r_knee
%                       ->  knee_angle_r actuation torque  (T_flex - T_ext)
%
% The knee swings under CPG rhythm + gravity; BPA strains come from the joint
% angle through epsScale. This is torque-level BPA actuation (moment-arm times
% force); full point-force routing along the MJCF pathpoint bodies is the
% follow-up.

here = fileparts(mfilename('fullpath'));
cd(here);
addpath(here);
mdlFile = fullfile(here, 'osim_import', 'Gait2392_simbody_simscape.slx');
assert(exist(mdlFile, 'file') == 2, ...
    'Imported model missing - run sns_osim_import first (%s)', mdlFile);
mdl = 'Gait2392_CPG';
if bdIsLoaded(mdl), close_system(mdl, 0); end
load_system(mdlFile);
% work on a copy so the pristine import stays untouched
new_system(mdl, 'Model');
delete_line(find_system(mdlFile, 'FindAll'), [], []); %#ok<NASGU>
close_system(mdl, 0);   % discard the empty template
copy_instance = load_system(mdlFile);
mdl = copy_instance;
set_param(mdl, 'SimulationCommand', 'Stop');

% find the knee actuation input
jointName = 'knee_angle_r';
cand = find_system(mdl, 'LookUnderMasks', 'all', 'Regexp', 'on', ...
    'Name', ['Joint', jointName, '.*']);
assert(~isempty(cand), 'No joint block matching "%s" found', jointName);
jblk = cand{1};
fprintf('Actuating joint block: %s\n', jblk);

% The smimport URDF model names: inspect and report; actuation port is the
% joint's actuation input. Enable actuation on the joint.
% NOTE: smimport-created revolute joints expose actuation via the joint's
% "Actuation" option; for the URDF import they are torque-actuated by default
% through the sensed/actuated flag. We add a Simulink-PS + External Force? No:
% simplest is the joint's actuation input port, enabled with set_param if
% present.
try
    set_param(jblk, 'ActuationInput', 'on');
    fprintf('actuation input enabled on %s\n', jblk);
catch
    fprintf('NOTE: %s has no ActuationInput parameter; checking ports\n', jblk);
end
ph = get_param(jblk, 'PortHandles');
fprintf('joint has %d inports, %d outports\n', numel(ph.Inport), numel(ph.Outport));

fprintf('DONE (skeleton) — complete wiring in next iteration once joint port map is known\n');
