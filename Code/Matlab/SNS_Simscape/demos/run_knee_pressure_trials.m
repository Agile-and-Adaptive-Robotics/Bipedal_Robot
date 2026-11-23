function run_knee_pressure_trials(outputFolder, family, protocol)
%RUN_KNEE_PRESSURE_TRIALS Simulate the saved knee with matched pulse controls.
% Quick protocol plots and checks six seconds around the default 3 s tap.
% The manuscript protocol preserves the original 30 s onset and 50 s run.
if nargin < 2, family='all'; end
if nargin < 3, protocol='manuscript'; end
here=fileparts(mfilename('fullpath'));
addpath(fileparts(here),here);
if ~isfolder(outputFolder),mkdir(outputFolder);end
if strcmp(family,'all')
    run_knee_pressure_trials(outputFolder,'nonspiking',protocol);
    run_knee_pressure_trials(outputFolder,'spiking',protocol);
    return
end
assert(ismember(family,{'nonspiking','spiking'}),'Unknown knee controller.');
if strcmp(family,'nonspiking'),mdl='KneeReflexDemo';else,mdl='KneeReflexDemo_Spiking';end
load_system(fullfile(here,[mdl '.slx']));
cleanup=onCleanup(@()close_system(mdl,0));
if strcmp(protocol,'quick')
    names={'control','pulse'};
    onset=3; ending=6;
elseif strcmp(protocol,'manuscript')
    names={'control','pulse','off_control','off_pulse','fine_control','fine_pulse','half_pulse','low_supply_control','low_supply_pulse'};
    if strcmp(family,'spiking')
        onset=3; ending=6;
    else
        onset=30; ending=50;   % same timing as the archived knee-reflex pulse
    end
else
    error('Unknown protocol');
end
p=knee_pressure_parameters();
signals={'th','thd','stim','V_MN_ext','V_MN_flex','u_ext','u_flex','A_ext','A_flex','P_ext','P_flex','L_ext','L_flex','F_ext','F_flex','T_ext','T_flex','Ia_ext_c','Ib_ext_c'};
extra={'IN_Ia_ext_spk','IN_Ib_ext_spk'};
if strcmp(family,'spiking'),signals=[signals,extra];end
for k=1:numel(names)
    tag=names{k}; gain=double(~startsWith(tag,'off_'));
    amp=2*double(endsWith(tag,'pulse')); supply=p.pressureSupply;
    if startsWith(tag,'half'),amp=amp/2;end
    if startsWith(tag,'low_supply'),supply=0.90*supply;end
    maxStep=0.001;if startsWith(tag,'fine'),maxStep=maxStep/2;end
    fprintf('Running %s %s, %s, %.1f s...\n',family,tag,protocol,ending);
    in=Simulink.SimulationInput(mdl);
    in=in.setVariable('feedbackGain',gain).setVariable('tapAmplitude',amp);
    in=in.setVariable('pressureSupply',supply).setVariable('tapOnset',onset);
    in=in.setModelParameter('StopTime',num2str(ending),'MaxStep',num2str(maxStep), ...
        'RelTol','1e-8','AbsTol','1e-10','ReturnWorkspaceOutputs','on');
    ticRun=tic; result=sim(in); wallSeconds=toc(ticRun);
    time=(0:0.001:ending).'; samples=zeros(numel(time),numel(signals)); spikeTimes=struct;
    for i=1:numel(signals)
        ts=result.get(['log_' signals{i}]);
        [knots,ix]=unique(ts.Time,'last');
        values=double(ts.Data(ix));
        if endsWith(signals{i},'_spk')
            active=double(ts.Data(:))>0.5;
            spikeTimes.(signals{i})=ts.Time(active & ~[false;active(1:end-1)]);
        end
        if strcmp(signals{i},'stim') || endsWith(signals{i},'_spk')
            method='previous';
        else
            method='linear';
        end
        samples(:,i)=interp1(knots,values,time,method);
        assert(all(isfinite(samples(:,i))),'Nonfinite %s in %s',signals{i},tag);
    end
    dest=fullfile(outputFolder,[family '_' tag '.mat']);
    save(dest,'time','samples','signals','mdl','family','tag','protocol','p','onset','ending','amp','gain','supply','maxStep','wallSeconds','spikeTimes','-v7');
    fprintf('Saved %s in %.1f s. Final knee angle %.8f deg. Baseline range %.5g deg.\n', ...
      dest,wallSeconds,rad2deg(samples(end,1)),rad2deg(range(samples(time>onset-0.5&time<onset,1))));
end
end
