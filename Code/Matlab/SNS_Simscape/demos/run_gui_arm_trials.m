function run_gui_arm_trials(outputFolder)
% Compare the two saved realizations of the pneumatic arm with reflexes on/off.
here=fileparts(mfilename('fullpath')); addpath(here,fileparts(here));
if ~isfolder(outputFolder),mkdir(outputFolder);end
signals={'th','mCup','F_bi','F_tri','A_bi','A_tri','V_MN_bi','V_MN_tri','P_bi'};
for tag={'source','gui'}
    label=tag{1};src=strcmp(label,'source');if src,mdl='BeerCupReflexDemo';else,mdl='ArmReflexGUI';end
    load_system(fullfile(here,[mdl '.slx']));
    syn=find_system(mdl,'SearchDepth',1,'Regexp','on','Name','^syn_');assert(numel(syn)==6);
    for enabled=[1 0]
        in=Simulink.SimulationInput(mdl);in=in.setModelParameter('ReturnWorkspaceOutputs','on');
        if src
            in=in.setVariable('kReflex',enabled);
        elseif ~enabled
            for i=1:numel(syn),in=in.setBlockParameter(syn{i},'gmax','0');end
        end
        tic;out=sim(in);wallSeconds=toc;time=(0:.005:10).'; samples=zeros(numel(time),numel(signals));
        for i=1:numel(signals)
            ts=out.get(['log_' signals{i}]);[knots,ix]=unique(ts.Time,'last');
            samples(:,i)=interp1(knots,ts.Data(ix),time);assert(all(isfinite(samples(:,i))));
        end
        save(fullfile(outputFolder,sprintf('arm_%s_%d.mat',label,enabled)),'time','samples','signals','label','enabled','wallSeconds','-v7');
        fprintf('Arm %s feedback %d: sag %.6f degrees; biceps %.6f, triceps %.6f N.\n',label,enabled,rad2deg(samples(end,1)),samples(end,3),samples(end,4));
    end
    close_system(mdl,0);
end
end
