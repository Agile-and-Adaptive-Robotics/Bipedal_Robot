function add_knee_activation(mdl, name, initial, position)
% Bounded neural excitation is supplied by the MN voltage-ramp output.
% Appendix A Thelen activation/deactivation equation (seconds); the output
% becomes a normalized pneumatic pressure command, not a muscle force gain.
block=[mdl '/' name];
add_block('simulink/Ports & Subsystems/Subsystem',block,'Position',position);
delete_line(block,'In1/1','Out1/1');
add_block('simulink/Math Operations/Sum',[block '/Difference'],'Inputs','+-','Position',[95 60 120 90]);
add_block('simulink/Math Operations/Gain',[block '/ScaleActivation'],'Gain','1.5','Position',[135 175 180 200]);
add_block('simulink/Math Operations/Bias',[block '/TimeScale'],'Bias','0.5','Position',[205 175 249 200]);
add_block('simulink/Math Operations/Gain',[block '/TurnOn'],'Gain','activationTime','Position',[271 120 315 145]);
add_block('simulink/Sources/Constant',[block '/OffConstant'],'Value','deactivationTime','Position',[200 238 250 259]);
add_block('simulink/Math Operations/Product',[block '/TurnOff'],'Inputs','*/','Position',[271 191 315 221]);
add_block('simulink/Signal Routing/Switch',[block '/SelectTime'],'Criteria','u2 > Threshold','Threshold','0','Position',[345 114 385 183]);
add_block('simulink/Math Operations/Product',[block '/Rate'],'Inputs','*/','Position',[409 65 444 105]);
add_block('simulink/Continuous/Integrator',[block '/Activation'],'InitialCondition',initial,'LimitOutput','on','LowerSaturationLimit','0','UpperSaturationLimit','1','Position',[469 64 501 98]);
set_param([block '/In1'],'Position',[25 61 55 77]);set_param([block '/Out1'],'Position',[550 70 580 85]);
lines={'In1/1','Difference/1';'Activation/1','Difference/2';'Activation/1','ScaleActivation/1'; ...
'ScaleActivation/1','TimeScale/1';'TimeScale/1','TurnOn/1';'OffConstant/1','TurnOff/1';'TimeScale/1','TurnOff/2'; ...
'TurnOn/1','SelectTime/1';'Difference/1','SelectTime/2';'TurnOff/1','SelectTime/3'; ...
'Difference/1','Rate/1';'SelectTime/1','Rate/2';'Rate/1','Activation/1';'Activation/1','Out1/1'};
for k=1:size(lines,1),add_line(block,lines{k,1},lines{k,2},'autorouting','on');end
end
