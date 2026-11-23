function sns_build_biological_knee(outputFolder)
% GUI-generated feedback circuit connected to two Gait2392 Thelen2003 muscles.
% One-DOF knee fixture (positive flexion). All forces, lengths, activation,
% activation dynamics, and motion are computed forward in native Simulink.
here=fileparts(mfilename('fullpath'));
if nargin<1,outputFolder=fullfile(tempdir,'sns_gui_biological_knee');end
[mdl,~]=reflex_gui_network('KneeReflexDemo','KneeReflexHuman',outputFolder);
load_system('SNS_Library');
j=jsondecode(fileread(fullfile(here,'gait2392_reflex_geometry.json')));
muscleNames={'vas_med_r','bifemsh_r'}; sides={'ext','flex'};
a0=.12; q0=pi/4; eq=j.opensim_equilibria;
ma=zeros(1,2); force0=zeros(1,2); fiber0=zeros(1,2); pars=cell(1,2);
for i=1:2
    x=j.muscles.(muscleNames{i});ma(i)=x.length_gradient_m_per_rad(46);
    k=find([eq.angle]==45 & [eq.activation]==a0 & strcmp({eq.muscle},muscleNames{i}),1);
    force0(i)=eq(k).force_N; fiber0(i)=eq(k).fiber_length_m;
    pars{i}=[x.Fmax x.l_opt x.l_slack deg2rad(x.alpha_deg) x.v_max x.FmaxTendonStrain ...
        x.FmaxMuscleStrain x.KshapePassive x.KshapeActive x.Af x.Flen x.fv_linear_extrap_threshold x.minimum_activation];
end
balancedLoad=sum(force0.*ma); % external flexion moment balancing tonic co-contraction
% At 45 degrees both Ia cells receive 3 nA. Ib is proportional to each
% equilibrium force. Choose tonic input to keep both motor commands at a0.
ia=max(0,min(1,(-52+3/.4+45)/10));
ib=max(0,min(1,(-52+min(10,force0.*[.020 .0222])/.4+45)/10));
v0=-45+10*a0;
desc=.5*(v0+52)+.05*(.18*ia*v0+(.20*ib+.15*ia)*(v0+72));
params=struct('activationTime',j.muscles.vas_med_r.activation_time_constant, ...
    'deactivationTime',j.muscles.vas_med_r.deactivation_time_constant, ...
    'sensoryPulse',2,'tapTorque',0,'eventTime',20,'pulseDuration',.03,'torqueDuration',.025);
fn=fieldnames(params); preload='';
for i=1:numel(fn),preload=[preload sprintf('%s=%.16g;\n',fn{i},params.(fn{i}))];end %#ok<AGROW>
set_param(mdl,'SolverType','Variable-step','Solver','ode15s','MaxStep','0.001','RelTol','1e-7','AbsTol','1e-9', ...
    'StopTime','35','PreLoadFcn',preload,'SignalLogging','off','UnconnectedInputMsg','none', ...
    'UnconnectedOutputMsg','none','Description','GUI-built Ia/Ib circuit with two compliant Gait2392 musculotendon blocks and a one-degree-of-freedom knee.');
evalin('base',preload);
B=@(name,lib,pos,varargin)add_block(lib,[mdl '/' name],'Position',pos,varargin{:});
W=@(a,b)add_line(mdl,a,b,'autorouting','on');
% Positive flexion: I*qdd = external torque - sum(F*dL/dq) - b*qd - K*(q-q0).
B('Moment','simulink/Math Operations/Sum',[523 729 548 791],'Inputs','++----');
B('Inertia','simulink/Math Operations/Gain',[599 724 657 758],'Gain','1/0.12');
B('Speed','simulink/Continuous/Integrator',[685 723 722 757]);
B('Angle','simulink/Continuous/Integrator',[754 723 797 757],'InitialCondition',num2str(q0,17));
B('Damping','simulink/Math Operations/Gain',[740 795 800 825],'Gain','0.6');
B('RestAngle','simulink/Math Operations/Bias',[803 800 853 830],'Bias',num2str(-q0,17));
B('ReturnSpring','simulink/Math Operations/Gain',[871 800 924 830],'Gain','0.5');
B('Counterweight','simulink/Sources/Constant',[420 722 465 751],'Value',num2str(balancedLoad,17));
W('Moment/1','Inertia/1');W('Inertia/1','Speed/1');W('Speed/1','Angle/1');
W('Speed/1','Damping/1');W('Angle/1','RestAngle/1');W('RestAngle/1','ReturnSpring/1');
W('Counterweight/1','Moment/1');W('Damping/1','Moment/5');W('ReturnSpring/1','Moment/6');
addPulse(mdl,'IaPulse','sensoryPulse','pulseDuration',[10 70]);
addPulse(mdl,'FlexionTorque','tapTorque','torqueDuration',[150 650]);
W('FlexionTorque/1','Moment/2');
records={'Angle','th';'Speed','thd';'IaPulse','stim';'FlexionTorque','tap'};
for i=1:2
    s=sides{i}; x=j.muscles.(muscleNames{i}); y=90+260*(i-1);
    n=@(base)[base '_' s]; out=@(block)[block '/1'];
    l=n('PathLength'); r=n('PathGradient'); bio=muscleNames{i};
    grid=mat2str(deg2rad(j.flexion_deg.'),17);
    B(l,'simulink/Lookup Tables/1-D Lookup Table',[1680 y+125 1770 y+157], ...
      'BreakpointsForDimension1',grid,'Table',mat2str(x.length_m.',17),'InterpMethod','Linear point-slope','ExtrapMethod','Clip');
    B(r,'simulink/Lookup Tables/1-D Lookup Table',[1680 y+165 1770 y+197], ...
      'BreakpointsForDimension1',grid,'Table',mat2str(x.length_gradient_m_per_rad.',17),'InterpMethod','Linear point-slope','ExtrapMethod','Clip');
    W('Angle/1',[l '/1']);W('Angle/1',[r '/1']);
    a=n('Activation'); add_knee_activation(mdl,a,num2str(a0,17),[1515 y 1605 y+59]);
    set_param([mdl '/' a '/Activation'],'LowerSaturationLimit',num2str(x.minimum_activation,17));
    W([n('MN') '/2'],[a '/1']);
    addMuscle(mdl,bio,pars{i},fiber0(i),[1820 y 1935 y+80],here);
    W(out(a),[bio '/1']);W(out(l),[bio '/2']);
    moment=n('MuscleMoment');B(moment,'simulink/Math Operations/Product',[1990 y+50 2025 y+84]);
    W([bio '/1'],[moment '/1']);W(out(r),[moment '/2']);W(out(moment),sprintf('Moment/%d',2+i));
    % Physical spindle input from normalized fiber stretch and fiber velocity.
    len=n('FiberExcursion'); vel=n('FiberVelocity'); base=n('SpindleRest');
    B(len,'simulink/Math Operations/Gain',[80 y+75 136 y+105],'Gain',num2str(1/x.l_opt,17));
    B(base,'simulink/Math Operations/Bias',[157 y+75 212 y+105],'Bias',num2str(.5-fiber0(i)/x.l_opt,17));
    B(vel,'simulink/Math Operations/Gain',[155 y+116 217 y+145],'Gain',num2str(1/x.l_opt,17));
    iaBlock=n('Spindle'); ibBlock=n('TendonOrgan');
    B(iaBlock,'SNS_Library/IaMuscleSpindle',[260 y+25 330 y+90],'Wl','6','Wv','3','Imax','10');
    B(ibBlock,'SNS_Library/IbGolgiTendon',[260 y+101 330 y+163],'Kf',num2str(.02+.0022*(i-1),17),'Imax','10');
    W([bio '/2'],[len '/1']);W(out(len),[base '/1']);W([bio '/3'],[vel '/1']);
    W(out(base),[iaBlock '/1']);W(out(vel),[iaBlock '/2']);W([bio '/1'],[ibBlock '/1']);
    if i==1
        B('IaInjected','simulink/Math Operations/Sum',[408 y 442 y+36],'Inputs','++');
        W(out(iaBlock),'IaInjected/1');W('IaPulse/1','IaInjected/2');W('IaInjected/1',[n('SN_Ia') '/1']);
    else
        W(out(iaBlock),[n('SN_Ia') '/1']);
    end
    W(out(ibBlock),[n('SN_Ib') '/1']);
    c=n('Tonic');B(c,'simulink/Sources/Constant',[1125 y-16 1202 y+11],'Value',num2str(desc(i),17));
    W(out(c),[n('MN') '/1']);
    records=[records;{n('MN'),n('V_MN');a,n('A');bio,n('F');iaBlock,n('Ia');ibBlock,n('Ib'); ...
        l,n('Lmt');[bio '/2'],n('Lfiber');[bio '/3'],n('Vfiber');moment,n('T')}]; %#ok<AGROW>
end
for i=1:size(records,1)
    b=records{i,1}; if ~contains(b,'/'),b=[b '/1'];end
    dest=['log_' records{i,2}]; y=48*i;
    B(dest,'simulink/Sinks/To Workspace',[2160 y 2240 y+25],'VariableName',dest,'SaveFormat','Timeseries');
    W(b,[dest '/1']);
end
save_system(mdl,fullfile(here,[mdl '.slx']));close_system(mdl,0);close_system('KneeReflexDemo',0);
fprintf('Created %s with elastic Gait2392 antagonists, %.8f Nm static counter-torque.\n',mdl,balancedLoad);
end

function addPulse(mdl,name,amplitude,duration,xy)
x=xy(1);y=xy(2);on=[name 'On'];off=[name 'Off'];
add_block('simulink/Sources/Step',[mdl '/' on],'Position',[x y x+27 y+20],'Time','eventTime','After',amplitude);
add_block('simulink/Sources/Step',[mdl '/' off],'Position',[x y+30 x+27 y+50],'Time',['eventTime+' duration],'After',['-' amplitude]);
add_block('simulink/Math Operations/Sum',[mdl '/' name],'Position',[x+62 y+8 x+87 y+40],'Inputs','++');
add_line(mdl,[on '/1'],[name '/1'],'autorouting','on');add_line(mdl,[off '/1'],[name '/2'],'autorouting','on');
end

function addMuscle(mdl,name,params,initial,position,here)
b=[mdl '/' name];add_block('built-in/Subsystem',b,'Position',position);
B=@(n,src,pos,varargin)add_block(src,[b '/' n],'Position',pos,varargin{:});
W=@(src,dst)add_line(b,src,dst,'autorouting','on');
B('Activation','simulink/Sources/In1',[25 35 55 50],'Port','1');
B('PathLength','simulink/Sources/In1',[25 70 55 85],'Port','2');
B('MuscleParameters','simulink/Sources/Constant',[50 165 148 193],'Value',mat2str(params,17));
B('FiberLength','simulink/Continuous/Integrator',[352 134 400 171],'InitialCondition',num2str(initial,17), ...
  'LimitOutput','on','LowerSaturationLimit',num2str(max(.01*params(2),params(2)*sin(params(4))/sin(acos(.1))),17),'UpperSaturationLimit','inf');
B('Thelen','simulink/User-Defined Functions/MATLAB Function',[177 34 288 125]);
chart=find(sfroot,'-isa','Stateflow.EMChart','Path',[b '/Thelen']);
chart.Script=fileread(fullfile(here,'thelen_reflex_fiber.m'));
B('Force','simulink/Sinks/Out1',[467 40 502 53],'Port','1');
B('Fiber','simulink/Sinks/Out1',[467 150 502 165],'Port','2');
B('FiberSpeed','simulink/Sinks/Out1',[467 101 502 115],'Port','3');
W('Activation/1','Thelen/1');W('PathLength/1','Thelen/2');W('FiberLength/1','Thelen/3');W('MuscleParameters/1','Thelen/4');
W('Thelen/1','Force/1');W('Thelen/2','FiberLength/1');W('FiberLength/1','Fiber/1');W('Thelen/2','FiberSpeed/1');
mask=Simulink.Mask.create(b);mask.Type='Gait2392 Thelen2003 biological muscle';mask.Display=sns_icon_muscle('BIO');
mask.Description='Elastic tendon, active and passive fibers, and constant-height pennation driven by Gait2392 parameters.';
end
