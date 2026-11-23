function run_gui_knee_trials(outputFolder, protocol)
% Run the saved GUI circuit with both neural and mechanical perturbations.
if nargin<2,protocol='full';end
here=fileparts(mfilename('fullpath'));addpath(here,fileparts(here));
if ~isfolder(outputFolder),mkdir(outputFolder);end
mdl='KneeReflexHuman';load_system(fullfile(here,[mdl '.slx']));
cleanup=onCleanup(@()close_system(mdl,0));
labels={'control','ia_pulse','torque_pulse','off_control','off_ia','off_torque','half_ia','fine_control','fine_ia'};
if strcmp(protocol,'quick'),labels=labels(1:3);end
syn=find_system(mdl,'SearchDepth',1,'Regexp','on','Name','^syn_');
assert(numel(syn)==6);
state=reflex_equilibrium_state(mdl);
logs=find_system(mdl,'SearchDepth',1,'BlockType','ToWorkspace');names=get_param(logs,'VariableName');
signals=cellfun(@(n)n(5:end),names,'UniformOutput',false);
for k=1:numel(labels)
    tag=labels{k}; off=startsWith(tag,'off_'); current=0; torque=0; onset=20; finalTime=35; maxStep=.001;
    if strcmp(protocol,'quick'),onset=2;finalTime=4;end
    if contains(tag,'ia'),current=2;end
    if startsWith(tag,'half'),current=current/2;end
    if contains(tag,'torque'),torque=1;end
    if startsWith(tag,'fine'),maxStep=.0005;end
    in=Simulink.SimulationInput(mdl);in=in.setInitialState(state);
    in=in.setVariable('sensoryPulse',current).setVariable('tapTorque',torque).setVariable('eventTime',onset);
    if off
        for i=1:numel(syn),in=in.setBlockParameter(syn{i},'gmax','0');end
        % Match the same 12% tonic motor activation, force, and joint posture
        % after disconnecting all synapses. Each uncoupled cell is at -43.8 mV.
        for s={'ext','flex'},in=in.setBlockParameter([mdl '/Tonic_' s{1}],'Value','4.1');end
    end
    in=in.setModelParameter('StopTime',num2str(finalTime),'MaxStep',num2str(maxStep),'ReturnWorkspaceOutputs','on');
    tic;out=sim(in);wallSeconds=toc;time=(0:.0005:finalTime).'; samples=zeros(numel(time),numel(signals));
    for i=1:numel(names)
        ts=out.get(names{i});[knots,ix]=unique(ts.Time,'last');values=double(ts.Data(ix));
        how='linear';if ismember(signals{i},{'stim','tap'}),how='previous';end
        samples(:,i)=interp1(knots,values,time,how);assert(all(isfinite(samples(:,i))),'Nonfinite: %s',signals{i});
    end
    save(fullfile(outputFolder,[tag '.mat']),'time','samples','signals','onset','tag','off','maxStep','wallSeconds','current','torque','-v7');
    i=find(strcmp(signals,'th'));fprintf('GUI human %s: before %.8f deg, range %.8f..%.8f, final %.8f (%.1f s).\n', ...
        tag,rad2deg(interp1(time,samples(:,i),onset-.001)),rad2deg(min(samples(time>=onset,i))), ...
        rad2deg(max(samples(time>=onset,i))),rad2deg(samples(end,i)),wallSeconds);
end
end
