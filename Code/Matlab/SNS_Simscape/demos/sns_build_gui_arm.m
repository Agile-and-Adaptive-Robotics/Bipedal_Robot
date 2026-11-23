function sns_build_gui_arm(outputFolder)
% Couple the original antagonistic 20 mm pneumatic arm to the GUI-built
% six-neuron circuit. Copy only external mechanical, sensory, and logging
% blocks from the original demonstration. Neural connections stay as built.
here=fileparts(mfilename('fullpath'));
if nargin<1,outputFolder=fullfile(tempdir,'sns_gui_loaded_arm');end
source='BeerCupReflexDemo';
[mdl,~]=reflex_gui_network(source,'ArmReflexGUI',outputFolder);
old=find_system(source,'SearchDepth',1,'Type','Block');
for i=1:numel(old)
    name=get_param(old{i},'Name'); ref=get_param(old{i},'ReferenceBlock');
    if ismember(ref,{'SNS_Library/NonSpikingSynapse','SNS_Library/NonSpikingNeuron'}),continue;end
    add_block(old{i},[mdl '/' name]);
end
for i=1:numel(old)
    dest=old{i}; name=get_param(dest,'Name'); ref=get_param(dest,'ReferenceBlock');
    if strcmp(ref,'SNS_Library/NonSpikingSynapse'),continue;end
    p=get_param(dest,'PortHandles');
    for j=1:numel(p.Inport)
        ln=get_param(p.Inport(j),'Line');if ln<0,continue;end
        sh=get_param(ln,'SrcPortHandle'); sb=get_param(sh,'Parent');
        if strcmp(get_param(sb,'ReferenceBlock'),'SNS_Library/NonSpikingSynapse'),continue;end
        src=sprintf('%s/%d',get_param(sb,'Name'),get_param(sh,'PortNumber'));
        add_line(mdl,src,sprintf('%s/%d',name,j),'autorouting','on');
    end
end
set_param(mdl,'SolverType','Variable-step','Solver',get_param(source,'Solver'), ...
    'PreLoadFcn',get_param(source,'PreLoadFcn'),'StopTime',get_param(source,'StopTime'), ...
    'RelTol',get_param(source,'RelTol'),'AbsTol',get_param(source,'AbsTol'), ...
    'MaxStep',get_param(source,'MaxStep'),'UnconnectedOutputMsg','none','UnconnectedInputMsg','none', ...
    'Description','GUI-built stretch, force, and reciprocal feedback from the biceps and triceps sides, coupled to a pneumatic elbow resisting an increasing cup load.');
save_system(mdl,fullfile(here,[mdl '.slx']));close_system(mdl,0);close_system(source,0);
end
