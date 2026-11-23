function test_robot_bridge(repo, outputDirectory)
% Compare a nonspiking motor-neuron command with a zero-input replay.
% repo is the Bipedal_Robot checkout. Both outputs are saved in outputDirectory.
if nargin < 2
    outputDirectory = tempname;
end
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end
originalDirectory = pwd;
repo = char(java.io.File(repo).getCanonicalPath());
here = char(java.io.File(outputDirectory).getCanonicalPath());
returnToCaller = onCleanup(@() cd(originalDirectory));
cd(here);
sns=fullfile(repo,'Code/Matlab/SNS_Simscape');
bs=fullfile(sns,'mujoco_bridge/blockset/mujoco-simulink-blockset-main');
addpath(sns,fullfile(bs,'blocks'),fullfile(bs,'tools'),bs);
setenv('PATH',[fullfile(bs,'blocks') ';' fullfile(bs,'lib','win64','mujoco','bin') ';' getenv('PATH')]);
src=fullfile(repo,'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_robot/gait2392_robot_cvt3.xml');
modelText=fileread(src);
meshDir=strrep(fileparts(src),'\','/');
modelText=strrep(modelText,'file="Geometry\', ['file="' meshDir '/Geometry/']);
DOM=xmlread(src);
world=DOM.getElementsByTagName('worldbody').item(0);
jointNodes=world.getElementsByTagName('joint'); jointNames=cell(1,jointNodes.getLength);
obs={'route_length','route_force'};
tags={'<actuatorpos name="route_length" actuator="vas_med_r"/>','<actuatorfrc name="route_force" actuator="vas_med_r"/>'};
for k=1:numel(jointNames)
    jointNames{k}=char(jointNodes.item(k-1).getAttribute('name'));
    obs{end+1}=sprintf('q_%03d',k);
    tags{end+1}=sprintf('<jointpos name="%s" joint="%s"/>',obs{end},jointNames{k});
end
assert(isempty(regexp(modelText,'<\s*sensor\b','once')));
modelText=strrep(modelText,'</mujoco>',['<sensor>' strjoin(tags,newline) '</sensor>' newline '</mujoco>']);
xmlPath=fullfile(here,'robotbody.xml');
fid=fopen(xmlPath,'w'); fwrite(fid,modelText,'char'); fclose(fid);
acts=DOM.getElementsByTagName('actuator').item(0).getElementsByTagName('general');
names=cell(1,acts.getLength);
for k=1:numel(names), names{k}=char(acts.item(k-1).getAttribute('name')); end
chan=find(strcmp(names,'vas_med_r')); assert(isscalar(chan));
uvec=zeros(numel(names),1); uvec(chan)=1; 
modelName='robotbody_sns_trial'; new_system(modelName); load_system('SNS_Library_R2025a'); modelWorkspace=get_param(modelName,'ModelWorkspace'); assignin(modelWorkspace,'actuatorMap',uvec);
add_block('SNS_Library_R2025a/NonSpikingNeuron',[modelName '/MotorNeuron'], ...
    'Vrest','0','Gm','1','Cm','10','Thr','0','Slope','5','Position',[100 50 215 115]);
add_block('simulink/Sources/Constant',[modelName '/InputCurrent'], 'Value','0','Position',[22 73 53 98]);
add_block('simulink/Math Operations/Gain',[modelName '/SelectMuscle'],'Gain','actuatorMap','Multiplication','Matrix(K*u)','Position',[290 86 355 120]);
add_block('simulink/Sinks/Terminator',[modelName '/Voltage'],'Position',[258 34 276 51]);
add_block('mjLib/MuJoCo Plant',[modelName '/RobotMechanics'],'xmlFileRel',strrep(xmlPath,'\','/'),'renderingType','None','rgbOutOption','off','depthOutOption','off','Orientation','left','Position',[175 192 352 250]);
add_block('simulink/Signal Routing/Bus Selector',[modelName '/SelectForce'],'OutputSignals','route_force','Orientation','left','Position',[113 205 128 238]);
add_block('simulink/Sinks/Out1',[modelName '/Readback'],'Orientation','left','Position',[28 211 56 231]);
add_line(modelName,'InputCurrent/1','MotorNeuron/1','autorouting','smart');
add_line(modelName,'MotorNeuron/1','Voltage/1','autorouting','smart');
add_line(modelName,'MotorNeuron/2','SelectMuscle/1','autorouting','smart');
add_line(modelName,'SelectMuscle/1','RobotMechanics/1','autorouting','smart');
add_line(modelName,'RobotMechanics/1','SelectForce/1','autorouting','smart');
add_line(modelName,'SelectForce/1','Readback/1','autorouting','smart');
ph=get_param([modelName '/MotorNeuron'],'PortHandles'); set_param(ph.Outport(2),'DataLogging','on','DataLoggingNameMode','Custom','DataLoggingName','motor');
ph=get_param([modelName '/RobotMechanics'],'PortHandles'); set_param(ph.Outport(1),'DataLogging','on','DataLoggingNameMode','Custom','DataLoggingName','readback');
set_param(modelName,'SolverType','Fixed-step','Solver','ode1','FixedStep','0.001','StopTime','0.050','SignalLogging','on','SignalLoggingName','trace');
result=struct('jointNames',{jointNames},'channel',chan,'allChannelCount',numel(names),'model','gait2392_robot_cvt3.xml','motor_neuron','SNS_Library_R2025a/NonSpikingNeuron','actuator','vas_med_r');
for k=1:2
    set_param([modelName '/InputCurrent'],'Value',num2str(2*(k-1)));
    res=sim(modelName,'ReturnWorkspaceOutputs','on'); bus=res.trace.get('readback').Values;
    cmd=res.trace.get('motor').Values;
    rec=struct(); rec.t=bus.route_length.Time(:); rec.command=cmd.Data(:);
    rec.force=bus.route_force.Data(:); rec.length=bus.route_length.Data(:); rec.q=zeros(numel(rec.t),numel(jointNames));
    for j=1:numel(jointNames), rec.q(:,j)=bus.(obs{j+2}).Data(:); end
    assert(all(isfinite([rec.force;rec.length;rec.q(:);rec.command])));
    fields={'off','on'}; result.(fields{k})=rec;
end

assert(max(result.on.command)>0.2 && max(abs(result.off.command))<1e-9);
result.delta_force=max(abs(result.on.force-result.off.force));
fid=fopen(fullfile(here,'bridge_record.json'),'w'); fwrite(fid,jsonencode(result),'char'); fclose(fid);
assert(result.delta_force>0.01,'Measured command effect did not exceed 0.01 N.');
fprintf('EXPERIMENT OK: %d actuators, channel %d, max delta force %.8g N, final S %.4g, samples %d\n',numel(names),chan,result.delta_force,result.on.command(end),numel(result.on.t));
set_param(modelName,'ScreenColor','white');
blocks=find_system(modelName,'SearchDepth',1,'Type','Block');
for b=1:numel(blocks), set_param(blocks{b},'FontName','Arial','FontSize','17','Selected','off'); end
set_param([modelName '/MotorNeuron'],'Name','SNS motor neuron');
set_param([modelName '/InputCurrent'],'Name','Current');
set_param([modelName '/SelectMuscle'],'Name','Actuator inputs');
set_param([modelName '/Voltage'],'ShowName','off');
set_param([modelName '/RobotMechanics'],'Name','MuJoCo mechanics');
set_param([modelName '/SelectForce'],'ShowName','off');
set_param([modelName '/Readback'],'Name','Force');
set_param(modelName,'SimulationCommand','Update');
hilite_system(modelName,'none');
% Use black ink for the library's yellow outline and name when printed.
set_param([modelName '/MuJoCo mechanics'],'ForegroundColor','black');
save_system(modelName,fullfile(here,[modelName '.slx']));
print(['-s' modelName],'-dpdf','-bestfit',fullfile(here,'robot_interface.pdf'));
print(['-s' modelName],'-dpng','-r300',fullfile(here,'robot_interface.png'));
close_system(modelName,0);
end
