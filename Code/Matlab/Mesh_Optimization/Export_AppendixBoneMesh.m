function Export_AppendixBoneMesh()
% Export_AppendixBoneMesh  Dissertation Appendix C static bone-mesh figure.
%
% Runs the STATIC bone-mesh section of Results/Knee_Extensor_20mm.m (the
% block ending in run('MuscleBonePlotting.m') -> AnimateKneeBoneMuscle
% 'Static',true) from a temporary copy of that script with only the GIF
% animation block (%% Moving muscle/bone plot ... clear bonePlotArgs)
% removed, then exports the returned static figure. The temporary copy is
% identical to the script of record otherwise: same dated result mat
% (Vas_Pam_20mm_Result_20260925.mat, set at its line 56), same context
% rebuild, same route multi-panel, same bonePlotArgs. The animation block
% is skipped so the existing Knee_Extensor_20mm.gif is never rewritten and
% no getframe/GIF capture runs headless; the script's local functions
% (below the animation block) are kept, because the static section calls
% them.
%
% Run headless:
%   matlab -batch "cd('D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization'); Export_AppendixBoneMesh"
%
% Outputs (Figures/94-AppendixC):
%   bone_mesh_static_extensor20mm.png  (200 dpi raster)
%   bone_mesh_static_extensor20mm.alt.txt (alt text sidecar)

%% Paths (recomputed after the run too: the analyzed script opens with
%% `clear`, which wipes this function workspace when run() returns)
thisDir = fileparts(mfilename('fullpath'));   % ...\Code\Matlab\Mesh_Optimization
resultsDir = fullfile(thisDir, 'Results');
srcFile = fullfile(resultsDir, 'Knee_Extensor_20mm.m');

%% Build the truncated script: drop ONLY the animation block
txt = fileread(srcFile);
lines = regexp(txt, '\r?\n', 'split');
lines = strrep(lines, sprintf('\r'), '');

iStatic = find(strcmp(strtrim(lines), 'run(''MuscleBonePlotting.m'')'), 1);
iClear = find(strcmp(strtrim(lines), 'clear bonePlotArgs'), 1);
assert(numel(iStatic) == 1 && numel(iClear) == 1, ...
    'Anchor lines not found uniquely in Knee_Extensor_20mm.m')
assert(iClear > iStatic, 'Unexpected block order in Knee_Extensor_20mm.m')
keep = [1:iStatic, iClear+1:numel(lines)];
fprintf('Truncating Knee_Extensor_20mm.m: dropping lines %d..%d (%d lines, animation block only).\n', ...
    iStatic+1, iClear, iClear-iStatic);

tempRun = fullfile(resultsDir, 'tmp_AppendixBoneMesh_run.m');
fid = fopen(tempRun, 'w');
fprintf(fid, '%s\n', lines{keep});
fclose(fid);

%% Run the static section (executes clc/clear/close all in THIS workspace,
%% wiping every local above; nothing below may rely on pre-run variables)
try
    run(tempRun)
catch err
    tempRun2 = fullfile(fileparts(mfilename('fullpath')), 'Results', 'tmp_AppendixBoneMesh_run.m');
    if isfile(tempRun2), delete(tempRun2); end
    rethrow(err)
end
tempRun2 = fullfile(fileparts(mfilename('fullpath')), 'Results', 'tmp_AppendixBoneMesh_run.m');
if isfile(tempRun2), delete(tempRun2); end

%% Locate the static figure (name set by AnimateKneeBoneMuscle) and export
figs = findobj(groot, 'Type', 'figure', ...
    'Name', 'Knee bones and muscle routes - static');
assert(~isempty(figs), 'Static bone-mesh figure was not created.')
fig = figs(end);

ax = gca(fig);
fprintf('Static bone figure: title "%s"\n', ax.Title.String);
fprintf('xBest p1 = [%.4f %.4f %.4f] m, pEnd = [%.4f %.4f %.4f] m, rest = %.4f m, tendon = %.4f m\n', ...
    p1, pEnd, rest, tendon);

root = fileparts(mfilename('fullpath'));   % ...\Code\Matlab\Mesh_Optimization (recomputed post-clear)
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', thisDir)
    end
    root = parent
end
figDir = fullfile(root, 'Documentation', 'Reports and Papers', ...
    'Dissertation', 'Figures', '94-AppendixC');
pngTarget = fullfile(figDir, 'bone_mesh_static_extensor20mm.png');
altTarget = fullfile(figDir, 'bone_mesh_static_extensor20mm.alt.txt');

exportgraphics(fig, pngTarget, 'Resolution', 200);

altText = [ ...
    'ALT TEXT (bone_mesh_static_extensor20mm.png): Static three-dimensional bone-mesh render of the full skeleton at exactly zero knee angle with the optimized 20 mm vasti extensor BPA route of the 2026-09-25 design of record (Vas_Pam_20mm_Result_20260925.mat). The axis title reads Knee angle = 0.0 deg, BPA path length = 0.784 m. Light-blue and near-black dot clouds are the OpenSim femur and tibia bone surfaces; blue clouds show the pelvis and hip, sacrum, lumbar spine, and foot. The magenta circled polyline is the optimized p1-to-pEnd BPA/tendon route with its origin (p1) and distal ring insertion (pEnd) endpoints circled; the light-orange polylines are the human Vastus Intermedius, Vastus Lateralis, and Vastus Medialis reference paths from OpenSim. Axes read X, m; Y, m; Z, m, and the legend lists all twelve entries.'];
fid = fopen(altTarget, 'w');
fprintf(fid, '%s', altText);
fclose(fid);

fprintf('Exported:\n  %s\n  %s\n', pngTarget, altTarget);
assert(isfile(pngTarget) && dir(pngTarget).bytes > 0, 'PNG export failed')
fprintf('Export_AppendixBoneMesh: done.\n')

end
