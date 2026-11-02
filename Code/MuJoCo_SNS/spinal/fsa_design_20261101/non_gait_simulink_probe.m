function non_gait_simulink_probe(libraryFile,recordDirectory)
% Isolated source SNS blocks only. No Simscape body or library regeneration.
here = fileparts(mfilename('fullpath'));
repo = fileparts(fileparts(fileparts(fileparts(here))));
if nargin < 1 || isempty(libraryFile)
    libraryFile = fullfile(repo,'Code','Matlab','SNS_Simscape','SNS_Library.slx');
end
if nargin < 2 || isempty(recordDirectory)
    recordDirectory = 'D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical';
end
if ~isfolder(recordDirectory), mkdir(recordDirectory); end
if bdIsLoaded('SNS_Library')
    loaded = get_param('SNS_Library','FileName');
    assert(strcmpi(loaded,libraryFile), ...
        'Different SNS_Library already loaded: %s. Preserve it and use a separate MATLAB session.',loaded);
end
record = struct('machine',getenv('COMPUTERNAME'),'matlab_release',version, ...
    'source_library',libraryFile,'physics_steps',0,'training_steps',0, ...
    'original_library_edited',false,'graded_cases',[],'spiking_cases',[]);
destination = fullfile(recordDirectory,'NON_GAIT_SIMULINK_BLOCK_CHECK.json');
assert(~isfile(destination), 'Preserve existing evidence: %s. Choose a new explicit recordDirectory.',destination);
originalLibraryBytes = readBytes(libraryFile);
try
    load_system(libraryFile);
    for tau = [0.03 0.015]
        mdl = sprintf('fsa_graded_port_%d',round(1e6*tau));
        if bdIsLoaded(mdl), error('Probe model name already in use: %s',mdl); end
        new_system(mdl);
        cleanup = onCleanup(@() closeOwn(mdl));
        configure(mdl,0.2);
        add_block('simulink/Sources/Constant',[mdl '/Current'],'Value','5');
        add_block('SNS_Library/NonSpikingNeuron',[mdl '/Relay'], ...
            'Vrest','0','Gm','1','Cm',num2str(1000*tau,17), ...
            'Thr','0','Slope','0.2');
        connect(mdl,'Relay');
        add_block('simulink/Sinks/To Workspace',[mdl '/RecordVm'], ...
            'VariableName','portVm','SaveFormat','Timeseries');
        add_line(mdl,'Relay/1','RecordVm/1');
        output = sim(mdl,'ReturnWorkspaceOutputs','on');
        trace = output.get('portVm');
        expected = 5*(1-exp(-trace.Time/tau));
        row = struct('desired_physical_tau_s',tau, ...
            'actual_mask_Cm_nF',str2double(get_param([mdl '/Relay'],'Cm')), ...
            'actual_mask_Gm_uS',str2double(get_param([mdl '/Relay'],'Gm')), ...
            'maximum_voltage_error_mV',max(abs(trace.Data(:)-expected(:))), ...
            'samples',numel(trace.Time),'final_voltage_mV',trace.Data(end));
        row.passed = row.maximum_voltage_error_mV < 1e-6;
        record.graded_cases = [record.graded_cases row];
        clear cleanup
    end
    % Current on-disk library may have been rebuilt without the additive
    % spiking blocks. Stage the EXACT source builder and a library copy.
    % Never unlock/regenerate the original library just to run this probe.
    if isempty(find_system('SNS_Library','SearchDepth',1,'Name','SpikingNeuron'))
        staged = fullfile(repo,'tmp','fsa_design_20261101_simulink');
        if ~isfolder(staged), mkdir(staged); end
        if bdIsLoaded('SNS_Library'), close_system('SNS_Library',0); end
        stagedLibrary = fullfile(staged,'SNS_Library.slx');
        stagedBuilder = fullfile(staged,'sns_build_spiking.m');
        copyfile(libraryFile,stagedLibrary);
        copyfile(fullfile(fileparts(libraryFile),'sns_build_spiking.m'),stagedBuilder);
        runStagedBuilder(stagedBuilder);
        record.spiking_probe_library = stagedLibrary;
        record.spiking_builder_copy = stagedBuilder;
        record.spiking_library_scope = 'Separate copied library plus original additive source builder; original installed library not regenerated.';
    else
        record.spiking_probe_library = libraryFile;
        record.spiking_library_scope = 'Explicit supplied library already contains source spiking blocks; no builder or library rewrite.';
    end
    for current = [2.9 5.3]
        mdl = sprintf('fsa_spike_port_%d',round(100*current));
        if bdIsLoaded(mdl), error('Probe model name already in use: %s',mdl); end
        new_system(mdl);
        cleanup = onCleanup(@() closeOwn(mdl));
        configure(mdl,0.6);
        add_block('simulink/Sources/Constant',[mdl '/Current'], ...
            'Value',num2str(current,17));
        add_block('SNS_Library/SpikingNeuron',[mdl '/Relay'], ...
            'Vrest','-52','Vth','-46','Vreset','-58','Gm','0.4','Cm','2');
        connect(mdl,'Relay');
        add_block('simulink/Sinks/To Workspace',[mdl '/RecordVm'], ...
            'VariableName','portVm','SaveFormat','Timeseries');
        add_line(mdl,'Relay/1','RecordVm/1');
        output = sim(mdl,'ReturnWorkspaceOutputs','on');
        trace = output.get('portVm');
        resets = find(diff(trace.Data(:)) < -6)+1;
        resetTimes = trace.Time(resets);
        resetTimes = resetTimes(resetTimes > 0.1);
        if numel(resetTimes) < 3, error('No repeated actual reset events recorded.'); end
        observedRate = 1/mean(diff(resetTimes));
        equilibrium = -52+current/0.4;
        analyticalRate = 1/(0.005*log((equilibrium+58)/(equilibrium+46)));
        row = struct('current_nA',current,'actual_mask_Cm_nF',2, ...
            'actual_mask_threshold_mV',str2double(get_param([mdl '/Relay'],'Vth')), ...
            'actual_mask_reset_mV',str2double(get_param([mdl '/Relay'],'Vreset')), ...
            'analytical_source_reset_rate_Hz',analyticalRate, ...
            'measured_reset_rate_Hz',observedRate, ...
            'relative_rate_error',abs(observedRate-analyticalRate)/analyticalRate, ...
            'reset_events_after_100ms',numel(resetTimes));
        row.passed = row.relative_rate_error < 0.02;
        record.spiking_cases = [record.spiking_cases row];
        clear cleanup
    end
    record.passed = all([record.graded_cases.passed]) && all([record.spiking_cases.passed]);
    record.original_library_bytes_unchanged = isequal(originalLibraryBytes,readBytes(libraryFile));
    record.passed = record.passed && record.original_library_bytes_unchanged;
    record.scope = ['Source physical-unit graded RC masks and fixed-threshold/' ...
        'reset-below-rest sensory rate ports. Not adaptive GLIF, whole walker, or Table2-complete family.'];
    writeRecord(destination,record);
    assert(record.passed,'An isolated electrical block check failed.');
catch exception
    record.passed = false;
    record.error = getReport(exception,'extended','hyperlinks','off');
    writeRecord(destination,record);
    rethrow(exception)
end
disp(destination);
end

function runStagedBuilder(path)
run(path);
end

function bytes = readBytes(path)
file = fopen(path,'r');
assert(file >= 0,'Cannot read source library.');
cleanup = onCleanup(@() fclose(file));
bytes = fread(file,Inf,'*uint8');
end

function configure(model,stop)
set_param(model,'StopTime',num2str(stop),'Solver','ode45', ...
    'MaxStep','1e-4','RelTol','1e-8','AbsTol','1e-9', ...
    'ReturnWorkspaceOutputs','on');
end

function connect(model,neuron)
add_line(model,'Current/1',[neuron '/1']);
add_block('simulink/Sources/Constant',[model '/ZeroSynapses'],'Value','[0;0]');
for port = 2:7
    add_line(model,'ZeroSynapses/1',sprintf('%s/%d',neuron,port),'autorouting','on');
end
end

function closeOwn(model)
if bdIsLoaded(model), close_system(model,0); end
end

function writeRecord(path,record)
file = fopen(path,'w','n','UTF-8');
assert(file >= 0,'Cannot write unit-check record.');
cleanup = onCleanup(@() fclose(file));
fwrite(file,jsonencode(record,'PrettyPrint',true),'char');
end
