function libraryFolder = prepare_spiking_library_extended(libraryFolder)
%PREPARE_SPIKING_LIBRARY_EXTENDED Add spiking blocks to a private library copy.
% Run in a fresh MATLAB batch process. The shared library is never saved.
here = fileparts(mfilename('fullpath'));
libraryFolder = char(libraryFolder);
if bdIsLoaded('SNS_Library')
    error('SNS:LoadedLibrary','Use a fresh MATLAB process; SNS_Library is already loaded.');
end
if isfile(fullfile(libraryFolder,'SNS_Library.slx'))
    error('SNS:ExistingLibrary','Choose a new workspace; its library will not be overwritten.');
end
if ~isfolder(libraryFolder), mkdir(libraryFolder); end
copyfile(fullfile(here,'SNS_Library.slx'),fullfile(libraryFolder,'SNS_Library.slx'));
copyfile(fullfile(here,'sns_build_spiking.m'),fullfile(libraryFolder,'sns_build_spiking.m'));
addpath(here);addpath(libraryFolder,'-begin');
prior = pwd; cleanup = onCleanup(@()cd(prior));
cd(libraryFolder);
runBuilder(fullfile(libraryFolder,'sns_build_spiking.m'));
load_system(fullfile(libraryFolder,'SNS_Library.slx'));
get_param('SNS_Library/SpikingNeuron','Handle');
get_param('SNS_Library/HybridSpikingSynapse','Handle');
fprintf('Private spiking library ready: %s\n',libraryFolder);
end

function runBuilder(path)
run(path);
end
