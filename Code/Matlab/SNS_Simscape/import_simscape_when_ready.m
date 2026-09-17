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

% Simscape Multibody Link export (Ben, 2026-09-16 17:59): 09_BA_003.xml —
% exports fine but DROPS the 4 SW Hinge mates (Hinge1/2/5/6, "not supported")
% and the 6 BPA mates to the assembly root ("components not resolved"); the
% PathMate drops are expected (patella paths). The knee DOF itself survives
% as primitive Concentric+Coincident mates between 04_02_KB_R_003 and
% 05_01_TI_R_006. Replace the SW Hinge mates with Concentric+Coincident and
% re-export to get the four link hinges as revolute joints.
% sw2urdf folder route (kept for reference):
%   ...\Knee assembly\09_BA_003.URDF\urdf\09_BA_003.URDF.urdf  (1-link skeleton)
inputPath = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Solid_Models\Biomimetics_2022-Knee_Test\Knee assembly\09_BA_003.xml';  % <-- edit
outDir = fileparts(mfilename('fullpath'));

% --- 0. license gate (either feature name counts) ----------------------------
licensed = (license('test','Simscape_Multibody') == 1) || (license('test','SimMechanics') == 1);
assert(licensed, ['Simscape Multibody is not licensed on this machine — but see the ' ...
    'LICENSE NOTE above; sns_urdf_smoke.m proves smimport runs under SimMechanics.']);

% --- 1. import URDF or Multibody-Link XML -> Simulink model ------------------
assert(isfile(inputPath), 'Input file not found: %s', inputPath);
% smimport derives the model name from the file name; '09_BA_003.URDF' has a
% dot (invalid) so MATLAB silently renames the model and smimport's return
% value is unusable. Import a sanitized COPY in the SAME folder — package://
% mesh paths resolve relative to the file, so the copy keeps its meshes.
srcDir = fileparts(inputPath);
[~, inExt, ~] = fileparts(inputPath);
tmpFile = fullfile(srcDir, ['mdl_knee_rig_import_tmp' inExt]);  % keep extension: smimport dispatches on it
if bdIsLoaded('mdl_knee_rig_import_tmp')
    close_system('mdl_knee_rig_import_tmp', 0);   % stale copy from a crashed run
end
copyfile(inputPath, tmpFile, 'f');
cleanup = onCleanup(@() delete(tmpFile));
% R2025b smimport's return value is NOT the model name (observed: a double),
% but with a sanitized file name the created model name is deterministic.
smimport(tmpFile);
model = 'mdl_knee_rig_import_tmp';
assert(bdIsLoaded(model), 'smimport did not create the expected model');
fprintf('Imported model: %s\n', model);

% --- 2. orientation / gravity ------------------------------------------------
% The CAD Z axis is up in 09_BA_003 (knee axis ~ Z, tibia pointing -Y).
set_param(model, 'StopTime', '2');
try
    gBlk = find_system(model, 'LookUnderMasks', 'all', ...
        'MaskType', 'Mechanism Configuration');
    if iscell(gBlk) && ~isempty(gBlk)
        % param is 'GravityVector' on R2025b (verified; 'Gravity' and
        % 'UniformGravity' are NOT settable — the latter is a mode dropdown)
        try
            cur = get_param(gBlk{1}, 'GravityVector');
        catch
            cur = '';
        end
        if ~strcmp(strtrim(cur), '[0 0 -9.80665]')
            set_param(gBlk{1}, 'GravityVector', '[0 0 -9.80665]');  % -Z if Z up in CAD
        end
        fprintf('gravity -Z confirmed on %s\n', gBlk{1});
    else
        fprintf('no Mechanism Configuration block found; gravity left at default\n');
    end
catch ME
    warning('gravity tweak failed: %s', ME.message);
end

% --- 3. save next to the SNS work --------------------------------------------
savePath = fullfile(outDir, [model '_imported.slx']);
save_system(model, savePath);
close_system(model, 0);
fprintf('Saved %s\n', savePath);

% --- 4. next steps (manual, in the GUI) --------------------------------------
% * URDF route: the sw2urdf joint list becomes revolute joints directly;
%   name the knee joint "knee" in the exporter dialog so it is easy to find.
% * XML route: Hinge mates (Hinge1/2/5/6) become revolute joints; the
%   "Knee Angle" angle mate carries the joint state target.
% * Add Joint Actuation (torque) to the knee revolute joint, driven by the
%   BPAForce outputs in KneeReflexDemo (replacing the reduced-order plant).
% * Check the Tibia / Theta1 Coordinate System frames for sensor frames.
