function outFile = run_sns_model_extended(modelName, outDir, stopTime, mode, libraryFolder)
%RUN_SNS_MODEL_EXTENDED Run an existing model without rebuilding or saving it.
% outFile = run_sns_model_extended('KneeReflexDemo',outDir,'5','normal')
% Solver and step size remain those of the saved model. Stop time and mode
% are temporary SimulationInput overrides. Results are never overwritten.
if nargin < 3, stopTime = '5'; end
if nargin < 4, mode = 'normal'; end
if nargin < 5, libraryFolder = ''; end
modelName = char(modelName); outDir = char(outDir);
if isnumeric(stopTime), stopTime = num2str(stopTime); end
if ~ismember(mode,{'normal','accelerator','rapid-accelerator'})
    error('SNS:Mode','Choose normal, accelerator, or rapid-accelerator.');
end
here = fileparts(mfilename('fullpath'));
addpath(here,fullfile(here,'demos'),fullfile(here,'results'));
if ~isempty(libraryFolder)
    addpath(libraryFolder,'-begin');
    load_system(fullfile(libraryFolder,'SNS_Library.slx'));
end
[~,mdl,ext] = fileparts(modelName);
if isempty(ext), mdl = modelName; modelPath = which([mdl '.slx']);
else, modelPath = modelName; end
if isempty(modelPath) || ~isfile(modelPath), error('SNS:MissingModel','Model file was not found.'); end
wasLoaded = bdIsLoaded(mdl);
if wasLoaded && strcmp(get_param(mdl,'Dirty'),'on')
    error('SNS:UnsavedModel','Use a saved model or a separate copy; this model has unsaved changes.');
end
if ~isfolder(outDir), mkdir(outDir); end
outFile = fullfile(outDir,[mdl '_run.mat']);
if isfile(outFile), error('SNS:ExistingRun','Choose a new output directory.'); end
load_system(modelPath);
cleanup = onCleanup(@()closeIfOpened(mdl,wasLoaded));
meta = struct('model',mdl,'modelPath',modelPath,'release',version, ...
    'solver',get_param(mdl,'Solver'),'fixedStep',get_param(mdl,'FixedStep'), ...
    'stopTime',char(stopTime),'mode',char(mode));
in = Simulink.SimulationInput(mdl);
in = in.setModelParameter('StopTime',char(stopTime),'SimulationMode',char(mode));
started = tic;
out = sim(in);
meta.wallSeconds = toc(started);
if bdIsLoaded('SNS_Library'), meta.libraryPath = get_param('SNS_Library','FileName'); end
save(outFile,'out','meta','-v7.3');
fprintf('Saved %s (%.3f s wall time)\n',outFile,meta.wallSeconds);
end

function closeIfOpened(mdl,wasLoaded)
if ~wasLoaded && bdIsLoaded(mdl), close_system(mdl,0); end
end
