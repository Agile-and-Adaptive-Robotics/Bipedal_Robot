function laptop_build_bridge()
% Build the MuJoCo-Simulink bridge natively on the laptop (R2025b + MinGW 8.1).
% Follows BRIDGE_REPORT.md "RESUMED & PROVEN" recipe: install.m with MJ_VER
% 2.3.7 -> tools\setupBuild (MINGW) -> mex -setup c++ -> tools\build.
% Assumes blockset\ was extracted from blockset.zip with:
%   - install.m MJ_VER edited to '2.3.7' (asserted here),
%   - tools\setupBuild.m selectedCompilerWin="MINGW" (asserted here),
%   - the AARL LOCAL PATCH applied to src\mj.cpp (asserted here) - upstream
%     initData() zero-fills qpos; without the keyframe reset our model
%     starts 0.95 m sunk in the floor and explodes.

here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', 'blockset', 'mujoco-simulink-blockset-main');

src = fileread(fullfile(repo, 'src', 'mj.cpp'));
assert(~isempty(regexp(src, 'AARL LOCAL PATCH', 'once')), ...
    'mj.cpp is missing the AARL LOCAL PATCH - apply it before building.');
inst = fileread(fullfile(repo, 'install.m'));
assert(~isempty(regexp(inst, 'MJ_VER = ''2.3.7''', 'once')), ...
    'install.m is not set to MJ_VER 2.3.7.');
sb = fileread(fullfile(repo, 'tools', 'setupBuild.m'));
assert(~isempty(regexp(sb, '^selectedCompilerWin="MINGW";', 'lineanchors', 'once')), ...
    'setupBuild.m is not switched to MINGW.');

cd(repo);
fprintf('--- install() (downloads mujoco 2.3.7 + glfw) ---\n');
install();
fprintf('--- setupBuild (MINGW) ---\n');
cd(fullfile(repo, 'tools'));
setupBuild();
fprintf('--- mex -setup c++ ---\n');
mex('-setup', 'c++');
fprintf('--- build (gmake) ---\n');
build;
blks = dir(fullfile(repo, 'blocks', '*.mexw64'));
names = {blks.name};
assert(~isempty(blks) && any(strcmp(names, 'mj_sfun.mexw64')) ...
    && any(strcmp(names, 'mj_sampletime.mexw64')), ...
    'gmake did not produce the expected mexw64 targets - see gmake output above.');
fprintf('=== BRIDGE BUILD COMPLETE (laptop R2025b): %s ===\n', strjoin(names, ', '));
end
