%% import_simscape_when_ready.m — one-command CAD import into Simscape Multibody
%
% CAD routes (either file works; set inputPath to whichever you produced):
%   A. SW2URDF exporter (PREFERRED, Ben 2026-09-09):
%      SolidWorks > Tools > Export > "Export as URDF..." (add-in sw2urdf 1.6.1;
%      installer already downloaded to Downloads\sw2urdfSetup_1.6.1.exe).
%      Save alongside the assembly, e.g.:
%      Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.URDF
%   B. Simscape Multibody Link (fallback; the add-in IS installed and
%      registered on this machine, just disabled at startup — enable in
%      SolidWorks Tools > Add-Ins): exports .xml + STEP parts.
%
% LICENSE NOTE (validated 2026-09-09 by sns_urdf_smoke.m): this license file
% reports license('test','Simscape_Multibody') = 0 but carries the product
% under the LEGACY feature name 'SimMechanics' (= 1). R2025b smimport runs
% fine — the old "Simscape Multibody not licensed" blocker is not real.
%
% Usage: set inputPath below, then run this script.

inputPath = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.URDF';  % <-- edit
outDir = fileparts(mfilename('fullpath'));

% --- 0. license gate (either feature name counts) ----------------------------
licensed = (license('test','Simscape_Multibody') == 1) || (license('test','SimMechanics') == 1);
assert(licensed, ['Simscape Multibody is not licensed on this machine — but see the ' ...
    'LICENSE NOTE above; sns_urdf_smoke.m proves smimport runs under SimMechanics.']);

% --- 1. import URDF or Multibody-Link XML -> Simulink model ------------------
assert(isfile(inputPath), 'Input file not found: %s', inputPath);
model = smimport(inputPath);
[~, modelBase, ~] = fileparts(inputPath);
fprintf('Imported model: %s\n', model);

% --- 2. orientation / gravity ------------------------------------------------
% The CAD Z axis is up in 09_BA_003 (knee axis ~ Z, tibia pointing -Y).
set_param(model, 'StopTime', '2');
try
    gBlk = find_system(model, 'LookUnderMasks', 'all', 'MaskType', 'Mechanism Configuration');
    if ~isempty(gBlk)
        set_param(gBlk{1}, 'Gravity', '[0 0 -9.80665]');   % -Z if Z up in CAD
        fprintf('gravity set to -Z on %s\n', gBlk{1});
    end
catch ME
    warning('gravity tweak failed: %s', ME.message);
end

% --- 3. save next to the SNS work -------------------------------------------
savePath = fullfile(outDir, [modelBase '_imported.slx']);
save_system(model, savePath);
fprintf('Saved %s\n', savePath);

% --- 4. next steps (manual, in the GUI) --------------------------------------
% * URDF route: the sw2urdf joint list becomes revolute joints directly;
%   name the knee joint "knee" in the exporter dialog so it is easy to find.
% * XML route: Hinge mates (Hinge1/2/5/6) become revolute joints; the
%   "Knee Angle" angle mate carries the joint state target.
% * Add Joint Actuation (torque) to the knee revolute joint, driven by the
%   BPAForce outputs in KneeReflexDemo (replacing the reduced-order plant).
% * Check the Tibia / Theta1 Coordinate System frames for sensor frames.
