function run_install()
% Run the blockset's install.m from its own folder (it uses relative paths).
% Tries savepath; if the default pathdef.m is not writable, falls back to
% userpath pathdef.m so the install still completes.
here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', 'blockset', 'mujoco-simulink-blockset-main');
cd(repo);
install();
end
