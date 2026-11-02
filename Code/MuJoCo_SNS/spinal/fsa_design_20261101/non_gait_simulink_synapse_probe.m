function non_gait_simulink_synapse_probe(tag,mode)
% Electrical event/chemical ports only; source library and models unchanged.
libraryFile = 'D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\simulink_neural_blocks\SNS_Library.slx';
if nargin<1,tag='';end
if nargin<2,mode='events';end
assert(all(isstrprop(tag,'alphanum') | tag=='_'),'Simple record tag required.');
destination = fullfile('D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical',['NON_GAIT_SYNAPSE_PORT_CHECK' tag '.json']);
assert(~isfile(destination),'Preserve previous evidence: %s',destination);
if bdIsLoaded('SNS_Library')
    assert(strcmpi(get_param('SNS_Library','FileName'),libraryFile),'Different source library already loaded.');
end
load_system(libraryFile);
record = struct('machine',getenv('COMPUTERNAME'),'source_library',libraryFile, ...
    'physics_steps',0,'training_steps',0,'source_primitive','additive/capped events, not reset-to-peak', ...
    'contract','100Hz events,10ms decay;mean signed current +/-1nA atpostVm-50mV', ...
    'cases',[]);
if strcmp(mode,'graded')
    record.source_primitive='graded source threshold/saturation conductance';
    record.contract='Independent signed endpoint current +/-1nA atpostVm-50mV;presynaptic domain-60..-40mV';
    record.cases=gradedPorts;
else
for hybrid = [false true]
    for reversal = [0 -80]
      for capped = [false true]
        model = sprintf('fsa_event_syn_%d_%d_%d',hybrid,abs(reversal),capped);
        assert(~bdIsLoaded(model),'Fixture already open.');new_system(model);
        cleanup = onCleanup(@() close_system(model,0));
        set_param(model,'StopTime','0.4','Solver','ode45','MaxStep','1e-4', ...
            'RelTol','1e-8','AbsTol','1e-10','ReturnWorkspaceOutputs','on');
        add_block('simulink/Sources/Pulse Generator',[model '/Events'], ...
            'Amplitude','1','Period','.01','PulseWidth','5','PhaseDelay','.02');
        source = 'Events/1';
        if hybrid
            add_block('simulink/Math Operations/Gain',[model '/VoltageScale'],'Gain','30');
            add_block('simulink/Math Operations/Bias',[model '/VoltageBase'],'Bias','-60');
            add_line(model,source,'VoltageScale/1');add_line(model,'VoltageScale/1','VoltageBase/1');source='VoltageBase/1';
        end
        block = 'SpikingSynapse';if hybrid, block='HybridSpikingSynapse';end
        increment = 1/(abs(reversal+50)*100*.01);
        if capped,increment=.2;end
        add_block(['SNS_Library/' block],[model '/Chemical'], ...
            'gmax','.12','ginc',num2str(increment,17),'tau_syn','10','Esyn',num2str(reversal));
        if hybrid, set_param([model '/Chemical'],'ThrPre','-45');end
        add_block('simulink/Sinks/To Workspace',[model '/Record'], ...
            'VariableName','chemicalPort','SaveFormat','Timeseries');
        add_line(model,source,'Chemical/1');add_line(model,'Chemical/1','Record/1');
        output=sim(model,'ReturnWorkspaceOutputs','on');trace=output.get('chemicalPort');
        t=trace.Time;values=squeeze(trace.Data);if size(values,1)~=numel(t),values=values';end
        expected=zeros(size(t));events=.02:.01:.4;previous=0;conductance=0;
        for event=events
            conductance=min(.12,conductance*exp(-(event-previous)/.01)+increment);
            selected=t>=event;expected(selected)=conductance*exp(-(t(selected)-event)/.01);previous=event;
        end
        away=true(size(t));for event=events,away=away & abs(t-event)>2e-6;end
        current=values(:,2)+50*values(:,1);tail=t>=.25;
        meanCurrent=trapz(t(tail),current(tail))/(t(find(tail,1,'last'))-t(find(tail,1)));
        target=1;if reversal<0,target=-1;end
        if capped,target=(reversal+50)*.12*100*.01*(1-exp(-1));end
        row=struct('hybrid',hybrid,'cap_regime',capped,'reversal_mV',reversal,'ginc_uS',increment, ...
            'tau_ms',10,'cap_uS',.12,'event_rate_Hz',100,'postVm_mV',-50, ...
            'max_conductance_error_uS',max(abs(values(away,1)-expected(away))), ...
            'signed_current_mean_nA',meanCurrent,'desired_mean_nA',target, ...
            'relative_mean_error',abs(meanCurrent-target)/abs(target), ...
            'mean_current_not_a_voltage_or_Fmax_claim',true);
        row.passed=row.max_conductance_error_uS<1e-5 && row.relative_mean_error<.02;
        record.cases=[record.cases row];clear cleanup;
      end
    end
end
end
record.passed=all([record.cases.passed]);
file=fopen(destination,'w');assert(file>=0);cleanup=onCleanup(@() fclose(file));
fwrite(file,jsonencode(record,PrettyPrint=true),'char');clear cleanup;
assert(record.passed,'Source event-port contract failed; diagnostic retained.');disp(destination);
end

function rows=gradedPorts
rows=[];times=(0:.01:.08)';voltage=[-70;-60;-55;-50;-45;-40;-30;-60;-60];
assignin('base','fsaPresynapticVoltage',[times voltage]);
for reversal=[0 -80]
    model=sprintf('fsa_graded_chemical_%d',abs(reversal));assert(~bdIsLoaded(model));new_system(model);
    cleanup=onCleanup(@() close_system(model,0));
    set_param(model,'StopTime','.075','Solver','ode45','MaxStep','1e-4','ReturnWorkspaceOutputs','on');
    add_block('simulink/Sources/From Workspace',[model '/Voltage'],'VariableName','fsaPresynapticVoltage','Interpolate','off','OutputAfterFinalValue','Holding final value');
    maximum=1/abs(reversal+50);
    add_block('SNS_Library/NonSpikingSynapse',[model '/Chemical'],'gmax',num2str(maximum,17),'Esyn',num2str(reversal),'ThrPre','-60','SlopePre','20');
    add_block('simulink/Signal Routing/Mux',[model '/Join'],'Inputs','2');
    add_block('simulink/Sinks/To Workspace',[model '/Record'],'VariableName','gradedPort','SaveFormat','Timeseries');
    add_line(model,'Voltage/1','Chemical/1');add_line(model,'Voltage/1','Join/1');add_line(model,'Chemical/1','Join/2');add_line(model,'Join/1','Record/1');
    output=sim(model,'ReturnWorkspaceOutputs','on');trace=output.get('gradedPort');values=squeeze(trace.Data);if size(values,1)~=numel(trace.Time),values=values';end
    fraction=min(1,max(0,(values(:,1)+60)/20));g=maximum*fraction;current=values(:,3)+50*values(:,2);wanted=fraction;if reversal<0,wanted=-wanted;end
    row=struct('reversal_mV',reversal,'gmax_uS',maximum,'ThrPre_mV',-60,'SlopePre_mV',20,'postVm_mV',-50,'fullscale_signed_current_nA',sign(reversal+50),'max_g_error_uS',max(abs(values(:,2)-g)),'maximum_signed_current_error_nA',max(abs(current-wanted)));
    row.passed=row.max_g_error_uS<1e-10 && row.maximum_signed_current_error_nA<1e-10;rows=[rows row];clear cleanup;
end
evalin('base','clear fsaPresynapticVoltage');
end
