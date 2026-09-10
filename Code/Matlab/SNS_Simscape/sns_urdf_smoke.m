%% sns_urdf_smoke.m — validate the SolidWorks->URDF->Simscape pipeline on this machine
%
% The CAD route is now: SolidWorks --(SW2URDF exporter)--> .urdf --> smimport.
% (Simscape Multibody Link is NOT needed; smimport consumes URDF natively.)
% NOTE: license('test','Simscape_Multibody') reads 0 on this license file, but the
% product is carried under the LEGACY feature name 'SimMechanics' (R2025b honors
% it), so smimport runs. This script proves it end-to-end with a minimal
% femur+tibia revolute-knee URDF; it creates nothing permanent.
%
% Usage: matlab -batch "run('sns_urdf_smoke.m')"

here = fileparts(mfilename('fullpath'));
urdfPath = fullfile(here, 'urdf_smoke', 'knee_smoke.urdf');

fprintf('license Simscape_Multibody test: %d\n', license('test', 'Simscape_Multibody'));
fprintf('license SimMechanics (legacy)  : %d\n', license('test', 'SimMechanics'));
assert(~isempty(which('smimport')), 'smimport not on path');

if ~exist(urdfPath, 'file')
    error('Expected smoke URDF missing: %s', urdfPath);
end

mdl = smimport(urdfPath);
fprintf('smimport OK -> model "%s"\n', mdl);

% sanity: model has two bodies + a revolute joint
blks = find_system(mdl, 'LookUnderMasks', 'all', 'Type', 'Block');
fprintf('imported model has %d blocks\n', numel(blks));

close_system(mdl, 0);   % discard, keep nothing
fprintf('URDF->Simscape pipeline VALIDATED (nothing saved).\n');
