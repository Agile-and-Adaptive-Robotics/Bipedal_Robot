%% sns_osim_import.m — import the MyoConverter Gait2392 MuJoCo model into Simscape
%
% Route (2026-09-10): OpenSim gait2392_simbody --MyoConverter--> MJCF (cvt3)
%   --mjcf2urdf.py--> URDF --smimport--> Simscape Multibody.
% The MJCF muscle pathpoints (slide-joint bodies) come along as links, so the
% later SNS/BPA wiring can attach to the same points the MuJoCo muscles use.
% Muscles/tendons themselves are NOT imported (URDF has no muscle element) —
% they become SNS_Library blocks driving the joints.
cdto = fileparts(mfilename('fullpath'));
cd(cdto);

urdf = fullfile(cdto, 'osim_import', 'gait2392_simbody.urdf');
assert(exist(urdf, 'file') == 2, 'URDF missing: %s', urdf);

fprintf('smimport on %s ...\n', urdf);
mdl = smimport(urdf);
fprintf('smimport OK -> model handle "%s"\n', mdl);

blks = find_system(mdl, 'LookUnderMasks', 'all', 'Type', 'Block');
fprintf('imported model has %d blocks\n', numel(blks));

% compile-check (updates the diagram without simulating)
set_param(mdl, 'StopTime', '0.01');
try
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('update (compile) OK\n');
    updateOK = true;
catch ME
    fprintf('update FAILED: %s\n', ME.message);
    updateOK = false;
end

outDir = fullfile(cdto, 'osim_import');
save_system(mdl, fullfile(outDir, 'Gait2392_simbody_simscape.slx'));
fprintf('saved %s\n', fullfile(outDir, 'Gait2392_simbody_simscape.slx'));
if updateOK
    fprintf('OSIM->SIMSCAPE IMPORT PASSED\n');
else
    fprintf('OSIM->SIMSCAPE IMPORT: imported + saved, but compile update failed (see above)\n');
end
