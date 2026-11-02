% Fresh verification of existing graded/spiking reflex models; source files unchanged.
here=fileparts(mfilename('fullpath')); repo=here;
while ~isfolder(fullfile(repo,'Code','Matlab','SNS_Simscape'))
    parent=fileparts(repo); assert(~strcmp(parent,repo)); repo=parent;
end
sns=fullfile(repo,'Code','Matlab','SNS_Simscape');
dest='D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\reflex_simulink';
if ~isfolder(dest), mkdir(dest); end
addpath(sns,fullfile(sns,'demos'));
Simulink.fileGenControl('set','CacheFolder',fullfile(dest,'cache'), ...
    'CodeGenFolder',fullfile(dest,'codegen'),'createDir',true);
disp(version); summary=struct();
% The current base library no longer contains the October spiking extension.
% Rebuild that additive extension in an isolated library snapshot, never the source.
libraryDir=fullfile(dest,'library'); if ~isfolder(libraryDir), mkdir(libraryDir); end
libraryCopy=fullfile(libraryDir,'SNS_Library.slx');
copyfile(fullfile(sns,'SNS_Library.slx'),libraryCopy);
load_system(libraryCopy);
assert(strcmp(get_param('SNS_Library','FileName'),libraryCopy));
run(fullfile(sns,'sns_build_spiking.m'));
summary.library_snapshot=libraryCopy;
for name={'KneeReflexDemo','KneeReflexDemo_Spiking'}
    mdl=name{1}; load_system(fullfile(sns,'demos',[mdl '.slx']));
    set_param(mdl,'StopTime','5'); timer=tic; out=sim(mdl); elapsed=toc(timer);
    theta=out.log_th; q=theta.Data(:)*180/pi; window=theta.Time>3;
    S=struct('model',mdl,'duration_s',theta.Time(end),'wall_including_compile_s',elapsed, ...
        'all_angles_finite',all(isfinite(q)),'knee_min_deg',min(q), ...
        'knee_max_deg',max(q),'final_deg',q(end),'settle_mean_deg',mean(q(window)));
    if strcmp(mdl,'KneeReflexDemo_Spiking')
        ia=out.log_IN_Ia_ext_spk; ib=out.log_IN_Ib_ext_spk;
        x=ia.Data(:)>=.5; y=ib.Data(:)>=.5;
        S.Ia_spikes=sum(x & ~[false;x(1:end-1)]);
        S.Ib_spikes=sum(y & ~[false;y(1:end-1)]);
        S.Ia_rate_Hz=S.Ia_spikes/5; S.Ib_rate_Hz=S.Ib_spikes/5;
    end
    assert(S.all_angles_finite && S.knee_min_deg>=0,'Knee failed finite/no-hyperextension gate');
    save(fullfile(dest,[mdl '.mat']),'out','S','-v7.3');
    summary.(mdl)=S; close_system(mdl,0);
    fid=fopen(fullfile(dest,'SUMMARY.json'),'w'); fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true)); fclose(fid);
    disp(S);
end
close_system('SNS_Library',0);
