function state=reflex_equilibrium_state(mdl)
% Set the six membrane voltages to their physiological tonic operating point
% through Simulink's state API, keeping the GUI-emitted blocks unchanged.
% Motor activation = 0.12; Ia/Ib resting inputs come from the physical model.
root=fileparts(mfilename('fullpath'));j=jsondecode(fileread(fullfile(root,'gait2392_reflex_geometry.json')));
eq=j.opensim_equilibria; q=eq([eq.angle]==45 & [eq.activation]==0.12);
ef=q(strcmp({q.muscle},'vas_med_r')).force_N;ff=q(strcmp({q.muscle},'bifemsh_r')).force_N;
volt=struct('MN_ext',-43.8,'MN_flex',-43.8,'SN_Ia_ext',-52+3/.4,'SN_Ia_flex',-52+3/.4, ...
    'SN_Ib_ext',-52+.02*ef/.4,'SN_Ib_flex',-52+.0222*ff/.4);
state=Simulink.BlockDiagram.getInitialState(mdl);
assert(isa(state,'Simulink.SimulationData.Dataset'),'Expected a native state dataset');
changed=0;
for i=1:state.numElements
    e=state.getElement(i);path=e.BlockPath.getBlock(1);names=strsplit(path,'/');
    if numel(names)==3 && strcmp(names{3},'membraneRC') && isfield(volt,names{2})
        e.Values.Data(:)=volt.(names{2});state=state.setElement(i,e);changed=changed+1;
    end
end
assert(changed==6,'Expected six initialized GUI membranes, found %d',changed);
end
