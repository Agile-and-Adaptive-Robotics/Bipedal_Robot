function run_transient(tag)
% Run from any directory; originals are never saved or rebuilt.
if nargin<1, tag='primary'; end
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here));
cd(here); diary(fullfile(here,[tag '_console.txt']));
disp(version); addpath(root); addpath(fullfile(root,'demos'));
load_system(fullfile(root,'demos','KneeReflexDemo.slx'));
evalin('base',get_param('KneeReflexDemo','PreLoadFcn'));
audit=dumpModel('KneeReflexDemo'); save('loaded_audit.mat','audit');
fid=fopen('loaded_audit.txt','w'); fprintf(fid,'%s',jsonencode(audit,PrettyPrint=true)); fclose(fid);
save_system('KneeReflexDemo',fullfile(here,'KneeTransient.slx'));
mdl='KneeTransient';
if ~strcmp(tag,'primary')
 for nm={'MN_ext','MN_flex'}, set_param([mdl '/' nm{1}],'Slope','10'); end
 for nm={'syn_Iaext_exc','syn_Iaflex_exc','syn_Iaext_inh_on_flex','syn_Iaflex_inh_on_ext'}
 set_param([mdl '/' nm{1}],'SlopePre','10');
 end
end
% Additional sensory current, downstream of spindle transducer, upstream of SN.
delete_line(mdl,'Ia_ext/1','SN_Ia_ext/1');
add_block('simulink/Math Operations/Sum',[mdl '/Ia_inject'],'Inputs','++');
add_block('simulink/Sources/Step',[mdl '/TapOn'],'Time','30','Before','0','After','tapAmplitude');
add_block('simulink/Sources/Step',[mdl '/TapOff'],'Time','30.03','Before','0','After','-tapAmplitude');
add_block('simulink/Math Operations/Sum',[mdl '/Tap'],'Inputs','++');
add_line(mdl,'TapOn/1','Tap/1'); add_line(mdl,'TapOff/1','Tap/2');
add_line(mdl,'Tap/1','Ia_inject/2'); add_line(mdl,'Ia_ext/1','Ia_inject/1');
add_line(mdl,'Ia_inject/1','SN_Ia_ext/1');
extra={'Tap',1,'stim';'SN_Ia_ext',1,'V_Ia_ext';'SN_Ia_flex',1,'V_Ia_flex';'SN_Ib_ext',1,'V_Ib_ext';'SN_Ib_flex',1,'V_Ib_flex';'MN_ext',2,'S_ext';'MN_flex',2,'S_flex';'syn_Iaext_exc',1,'syn_ext';'syn_Iaext_inh_on_flex',1,'syn_flex'};
for k=1:size(extra,1)
 add_block('simulink/Sinks/To Workspace',[mdl '/log_' extra{k,3}], 'VariableName',['log_' extra{k,3}], 'SaveFormat','Timeseries');
 add_line(mdl,sprintf('%s/%d',extra{k,1},extra{k,2}),['log_' extra{k,3} '/1']);
end
set_param(mdl,'Solver','ode45','RelTol','1e-8','AbsTol','1e-10','MaxStep','0.001','StopTime','50','ReturnWorkspaceOutputs','on');
set_param(mdl,'PreLoadFcn',[get_param(mdl,'PreLoadFcn') newline 'tapAmplitude=2;']);
assignin('base','tapAmplitude',0); save_system(mdl);
control=sim(mdl); save([tag '_control.mat'],'control','-v7.3');
reportBaseline(control);
if ~strcmp(tag,'primary')
 s=control.log_th; v=control.log_thd; ix=s.Time>=25&s.Time<30;
 assert(range(s.Data(ix))<1e-5 && max(abs(v.Data(ix)))<1e-5,'Operating point not settled');
 for nm={'A_ext','A_flex','V_MN_ext','V_MN_flex','V_Ia_ext','V_Ia_flex','V_Ib_ext','V_Ib_flex'}
 s=control.get(['log_' nm{1}]); ix=s.Time>=25&s.Time<30; assert(range(s.Data(ix))<1e-5,'Neural/activation state not settled');
 end
end
assignin('base','tapAmplitude',2); pulse=sim(mdl); save([tag '_pulse.mat'],'pulse','-v7.3');
set_param(mdl,'MaxStep','0.0005'); fine=sim(mdl); save([tag '_fine.mat'],'fine','-v7.3');
parameters=dumpModel(mdl); save([tag '_parameters.mat'],'parameters','extra');
quantify(control,pulse,fine,tag);
close_system(mdl,0); diary off;
end
function a=dumpModel(m)
b=find_system(m,'LookUnderMasks','all','FollowLinks','on','Type','Block'); a=struct;
for k=1:numel(b)
 p=get_param(b{k},'DialogParameters'); q=struct; if ~isempty(p), n=fieldnames(p); for j=1:numel(n), try,q.(n{j})=get_param(b{k},n{j});catch,end,end,end
 a.blocks(k).path=b{k}; a.blocks(k).parameters=q;
 ph=get_param(b{k},'PortConnectivity'); a.blocks(k).connectivity=ph;
end
a.preload=get_param(m,'PreLoadFcn');
end
function reportBaseline(c)
names={'th','thd','A_ext','A_flex','V_MN_ext','V_MN_flex','V_Ia_ext','V_Ia_flex','V_Ib_ext','V_Ib_flex','S_ext','S_flex'};
for k=1:numel(names)
 s=c.get(['log_' names{k}]); ix=s.Time>=25 & s.Time<30;
 fprintf('BASE %s mean %.12g range %.12g\n',names{k},mean(s.Data(ix)),range(s.Data(ix)));
end
end
function quantify(c,p,f,tag)
t=(0:0.001:50)'; names={'th','thd','V_MN_ext','V_MN_flex','A_ext','A_flex','F_ext','F_flex','stim'};
for k=1:numel(names)
 cp=c.get(['log_' names{k}]); pp=p.get(['log_' names{k}]); fp=f.get(['log_' names{k}]);
 C(:,k)=interp1(cp.Time,cp.Data,t); P(:,k)=interp1(pp.Time,pp.Data,t); F(:,k)=interp1(fp.Time,fp.Data,t);
end
ix=t>=30; base=mean(C(t>=25&t<30,1)); d=(P(:,1)-C(:,1))*180/pi;
[v,j]=min(d(ix)); tt=t(ix); summary.extension_deg=v; summary.peak_latency_s=tt(j)-30;
summary.flexion_peak_deg=max(d(ix)); summary.baseline_deg=base*180/pi;
summary.residual_deg=mean(P(t>=45,1)-base)*180/pi;
summary.control_drift_deg=range(C(t>=25,1))*180/pi;
summary.step_error_deg=max(abs(P(:,1)-F(:,1)))*180/pi;
summary.recovered=abs(summary.residual_deg)<0.01 && max(abs(P(t>=45,2)))<1e-4;
summary.paper_baseline_deg=-summary.baseline_deg;
summary.paper_extension_excursion_deg=-summary.extension_deg;
summary.paper_residual_deg=-summary.residual_deg;
bad=find(t>=30.03 & (abs((P(:,1)-base)*180/pi)>0.01 | abs(P(:,2))>1e-4),1,'last');
summary.recovery_s=t(bad+1)-30;
first=find(t>=30 & abs(d)>0.01,1); summary.onset_001deg_s=t(first)-30;
disp(summary); save([tag '_summary.mat'],'summary','t','C','P','F','names');
fid=fopen([tag '_summary.json'],'w'); fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true)); fclose(fid);
fig=figure('Visible','off','Position',[100 100 1000 1100]); tiledlayout(6,1);
sel=t>=29&t<=40; x=t(sel)-30;
ys={P(sel,9),-[P(sel,1),C(sel,1)]*180/pi,-P(sel,2),P(sel,3:4),P(sel,5:6),P(sel,7:8)};
labs={'Ia pulse (nA)','Knee angle (deg; flexion -)','Velocity (rad/s; flexion -)','MN voltage (mV)','Activation (0-1)','BPA force (N)'};
for k=1:6
 nexttile; plot(x,ys{k}(:,1),'-','LineWidth',1.2); hold on;
 if size(ys{k},2)>1, plot(x,ys{k}(:,2),'--','LineWidth',1.2); end
 ylabel(labs{k}); grid on; xline(0,':'); xlim([-1 10]);
 if k==2,legend('Stimulus','No stimulus');elseif k>=4,legend('Extensor','Flexor');end
end
xlabel('Time from pulse onset (s)'); set(findall(fig,'-property','FontName'),'FontName','Arial'); set(findall(fig,'-property','FontSize'),'FontSize',11);
exportgraphics(fig,[tag '_response.pdf'],'ContentType','vector'); exportgraphics(fig,[tag '_response.png'],'Resolution',180); savefig(fig,[tag '_response.fig']); close(fig);
end
