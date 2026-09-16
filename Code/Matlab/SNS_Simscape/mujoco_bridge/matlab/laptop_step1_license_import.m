function laptop_step1_license_import()
% Laptop (DESKTOP-5Q16KE9, R2025b) bring-up, step 1 — three checks in one
% MATLAB startup:
%  (1) E0 rerun: does a HAND-BUILT Simscape Multibody model (primitive
%      blocks) build AND run here? EB475WS4 fails at block-ADD time; Ben's
%      two-cylinder-elbow reduced-order plant needs this machine to PASS.
%  (2) MEX C++ compiler configuration (needed before the bridge mex build).
%  (3) smimport of the REAL 09_BA_003 URDF export (currently a 1-link
%      skeleton — validates the pipe end-to-end + mesh resolution).

here = fileparts(mfilename('fullpath'));
addpath(here);

fprintf('=== (1) E0 hand-built Simscape Multibody ===\n');
e0_multibody_license();

fprintf('=== (2) MEX compiler configuration ===\n');
cc = mex.getCompilerConfigurations('C++', 'Selected');
if isempty(cc)
    fprintf('MEX C++: NOT CONFIGURED (MinGW-w64 support package needed)\n');
else
    fprintf('MEX C++ selected: %s %s\n', cc(1).Name, cc(1).Version);
end
fprintf('MW_MINGW64_LOC: "%s"\n', getenv('MW_MINGW64_LOC'));

fprintf('=== (3) smimport of the real knee-assembly URDF ===\n');
snsDir = fullfile(here, '..', '..');
run(fullfile(snsDir, 'import_simscape_when_ready.m'));
fprintf('=== step 1 done ===\n');
end
