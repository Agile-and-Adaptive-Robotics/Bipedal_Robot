function [mdl, J] = reflex_gui_network(template, mdl, outputFolder)
% Build the working source circuit via the palette's actual Build callback.
% The generated neural blocks and connections remain intact when mechanics
% are connected. The optional output directory avoids overwriting a user's
% last circuit in the GUI.
here=fileparts(mfilename('fullpath')); addpath(here,fileparts(here),fullfile(fileparts(here),'gui'));
if ~isfolder(outputFolder),mkdir(outputFolder);end
load_system(fullfile(here,[template '.slx']));
blocks=find_system(template,'SearchDepth',1,'Type','Block');
neurons=struct('name',{},'Vrest_mV',{},'Gm_uS',{},'Cm_nF',{},'Thr_mV',{},'Slope_mV',{});
synapses=struct('src',{},'dst',{},'g_uS',{},'Esyn_mV',{},'ThrPre_mV',{},'SlopePre_mV',{});
for i=1:numel(blocks)
    b=blocks{i}; ref=get_param(b,'ReferenceBlock');
    if strcmp(ref,'SNS_Library/NonSpikingNeuron')
        values=propertiesOf(b,{'Vrest','Gm','Cm','Thr','Slope'});
        neurons(end+1)=struct('name',get_param(b,'Name'),'Vrest_mV',values(1), ...
            'Gm_uS',values(2),'Cm_nF',values(3),'Thr_mV',values(4),'Slope_mV',values(5)); %#ok<AGROW>
    elseif strcmp(ref,'SNS_Library/NonSpikingSynapse')
        ports=get_param(b,'PortHandles');
        src=get_param(get_param(ports.Inport(1),'Line'),'SrcBlockHandle');
        dst=get_param(get_param(ports.Outport(1),'Line'),'DstBlockHandle');
        assert(isscalar(dst));
        values=propertiesOf(b,{'gmax','Esyn','ThrPre','SlopePre'});
        synapses(end+1)=struct('src',get_param(src,'Name'),'dst',get_param(dst,'Name'), ...
            'g_uS',values(1),'Esyn_mV',values(2),'ThrPre_mV',values(3),'SlopePre_mV',values(4)); %#ok<AGROW>
    end
end
assert(numel(neurons)==6 && numel(synapses)==6,'Unexpected source topology');
J=struct('neurons',neurons,'synapses',synapses);
source=fullfile(outputFolder,'reflex_source.json');
fid=fopen(source,'w');fwrite(fid,jsonencode(J,'PrettyPrint',true));fclose(fid);
app=SNSGui(); app.f.Visible='off'; closing=onCleanup(@()delete(app.f));
app.loadJson(source); app.cbBuild([],[],outputFolder);
emitted=jsondecode(fileread(fullfile(outputFolder,'circuit_last.json')));
for i=1:numel(J.neurons)
    assert(isequal(J.neurons(i),emitted.neurons(i)),'GUI changed a neuron property');
end
for i=1:numel(J.synapses)
    assert(isequal(J.synapses(i),emitted.synapses(i)),'GUI changed a synapse property');
end
load_system(fullfile(outputFolder,'SNS_GuiCircuit.slx'));
if bdIsLoaded(mdl),close_system(mdl,0);end
save_system('SNS_GuiCircuit',fullfile(here,[mdl '.slx']));
load_system(fullfile(here,[mdl '.slx']));
for i=1:numel(neurons)
    n=neurons(i); b=[mdl '/' n.name]; values=propertiesOf(b,{'Vrest','Gm','Cm','Thr','Slope'});
    assert(max(abs(values-[n.Vrest_mV n.Gm_uS n.Cm_nF n.Thr_mV n.Slope_mV]))<1e-10);
end
for i=1:numel(synapses)
    s=synapses(i); b=sprintf('%s/syn_%d_%s_to_%s',mdl,i,s.src,s.dst);
    pt=get_param(b,'PortHandles');
    assert(strcmp(get_param(get_param(get_param(pt.Inport(1),'Line'),'SrcBlockHandle'),'Name'),s.src));
    assert(strcmp(get_param(get_param(get_param(pt.Outport(1),'Line'),'DstBlockHandle'),'Name'),s.dst));
    assert(max(abs(propertiesOf(b,{'gmax','Esyn','ThrPre','SlopePre'})-[s.g_uS s.Esyn_mV s.ThrPre_mV s.SlopePre_mV]))<1e-10);
end
fprintf('GUI Build: six neuron masks, six signed connections, and every exported weight verified in %s.\n',mdl);
end

function values=propertiesOf(block,names)
values=zeros(size(names));
for i=1:numel(names)
    values(i)=slResolve(get_param(block,names{i}),block);
end
end
