% Explicit Ib-afferent current-clamp verification on existing reflex mechanics.
here=fileparts(mfilename('fullpath')); repo=here;
while ~isfolder(fullfile(repo,'Code','Matlab','SNS_Simscape')), repo=fileparts(repo); end
sns=fullfile(repo,'Code','Matlab','SNS_Simscape');
dest='D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\reflex_simulink';
addpath(sns,fullfile(sns,'demos')); addpath(fullfile(dest,'library'),'-begin');
Simulink.fileGenControl('set','CacheFolder',fullfile(dest,'cache'),'CodeGenFolder',fullfile(dest,'codegen'),'createDir',true);
load_system(fullfile(dest,'library','SNS_Library.slx'));
mdl='KneeReflexDemo_Spiking'; load_system(fullfile(sns,'demos',[mdl '.slx']));
delete_line(mdl,'Ib_ext/1','IN_Ib_ext/1');
add_block('simulink/Math Operations/Sum',[mdl '/IbCurrentSum'],'Inputs','+++');
add_block('simulink/Sources/Step',[mdl '/IbPulseOn'],'Time','2','Before','0','After','3','OutDataTypeStr','double');
add_block('simulink/Sources/Step',[mdl '/IbPulseOff'],'Time','2.1','Before','0','After','-3','OutDataTypeStr','double');
add_line(mdl,'Ib_ext/1','IbCurrentSum/1');add_line(mdl,'IbPulseOn/1','IbCurrentSum/2');add_line(mdl,'IbPulseOff/1','IbCurrentSum/3');add_line(mdl,'IbCurrentSum/1','IN_Ib_ext/1');
add_block('simulink/Sinks/To Workspace',[mdl '/IbCurrentProbe'],'VariableName','Ib_input_nA','SaveFormat','Timeseries');add_line(mdl,'IbCurrentSum/1','IbCurrentProbe/1');
set_param(mdl,'StopTime','5'); timer=tic; out=sim(mdl);elapsed=toc(timer);
q=out.log_th.Data(:)*180/pi;sp=out.log_IN_Ib_ext_spk;x=sp.Data(:)>=.5;events=x & ~[false;x(1:end-1)];
S=struct('scope','Ib current clamp +3nA at2.0–2.1s added to actual sensor current; no imposed mechanical load or source-file edit.', ...
 'actual_total_Ib_spikes',sum(events),'Ib_spikes_during_clamp',sum(events & sp.Time>=2 & sp.Time<2.1), ...
 'knee_min_deg',min(q),'knee_max_deg',max(q),'all_finite',all(isfinite(q)),'wall_including_compile_s',elapsed);
assert(S.Ib_spikes_during_clamp>0 && S.knee_min_deg>=0 && S.all_finite);
save(fullfile(dest,'Ib_current_clamp.mat'),'out','S','-v7.3');fid=fopen(fullfile(dest,'Ib_current_clamp.json'),'w');fprintf(fid,'%s',jsonencode(S,PrettyPrint=true));fclose(fid);disp(S);
close_system(mdl,0);close_system('SNS_Library',0);
