%% import_simscape_when_ready.m — run this once Simscape Multibody + the Link plugin exist
%
% Ben: this automates the Simscape Multibody import of 09_BA_003 as soon as two
% blockers are cleared (see README_SNS_Simscape.md):
%   1. MATLAB license includes Simscape Multibody  (license('test','Simscape_Multibody') == 1)
%   2. The "Simscape Multibody Link" SolidWorks add-in is installed and the assembly
%      was exported to XML (Simscape Multibody Link > Export in SolidWorks).
%
% Usage:  edit smXmlPath below, then run this script.

smXmlPath = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.xml';  % <-- export output
outDir = fileparts(mfilename('fullpath'));

% --- 0. license gate ---------------------------------------------------------
% This license file carries Simscape Multibody under its LEGACY feature name
% "SimMechanics" (product renamed in R2012; R2025b honors the old name).
licensed = (license('test','Simscape_Multibody') == 1) || (license('test','SimMechanics') == 1);
if ~licensed
    error(['Simscape Multibody is not licensed on this machine. ' ...
        'Fix the license first (see README_SNS_Simscape.md).']);
end

% --- 1. import CAD XML -> Simulink model -------------------------------------
model = smimport(smXmlPath);                    % raw import
% variable names for the generated model
[~, modelBase, ~] = fileparts(smXmlPath);
fprintf('Imported model: %s\n', model);

% --- 2. orientation / gravity -------------------------------------------------
% The CAD Z axis is up in 09_BA_003 (knee axis ~ Z, tibia pointing -Y).
% smimport keeps CAD axes; set gravity accordingly and fix the ground body.
set_param(model, 'StopTime', '2');
try
    % gravity block name varies; find it
    gBlk = find_system(model, 'LookUnderMasks', 'all', 'MaskType', 'Mechanism Configuration');
    if ~isempty(gBlk)
        set_param(gBlk{1}, 'Gravity', '[0 0 -9.80665]');   % -Z if Z up in CAD
        fprintf('gravity set to -Z on %s\n', gBlk{1});
    end
catch ME
    warning('gravity tweak failed: %s', ME.message);
end

% --- 3. save -----------------------------------------------------------------
save_system(model, fullfile(outDir, [modelBase '_imported.slx']));
fprintf('Saved %s\n', fullfile(outDir, [modelBase '_imported.slx']));

% --- 4. next steps (manual, in the GUI) --------------------------------------
% * smimport creates revolute joints from the Hinge mates (Hinge1/2/5/6) and a
%   "Knee Angle" state target from the angle-dimension mate.
% * Add Joint Actuation (torque) to the knee revolute joint, driven by the
%   BPAForce outputs in KneeReflexDemo (replace the reduced-order plant).
% * Check the Tibia Coordinate System / Theta1 Coordinate System frames that
%   already exist in the assembly for sensor frames.
