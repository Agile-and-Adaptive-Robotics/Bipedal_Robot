function verify_outputs
here=fileparts(mfilename('fullpath')); cd(here); diary('verification.txt');
disp(checkcode('run_transient.m','-id')); disp(checkcode('audit_builder.m','-id'));
a=load('graded_summary.mat'); b=load('replay_summary.mat');
fprintf('Fresh process full resampled pulse equality: %d; max angle error %.12g rad\n',isequaln(a.P,b.P),max(abs(a.P(:,1)-b.P(:,1))));
disp(b.summary);
c=load('replay_control.mat'); c=c.control;
fig=figure('Visible','off','Position',[100 100 1000 900]); tiledlayout(4,1);
nexttile; plot(c.log_th.Time,-c.log_th.Data*180/pi); ylabel('Angle (deg)'); title('Full no-stimulus settling record; flexion negative');
nexttile; plot(c.log_thd.Time,-c.log_thd.Data); ylabel('Velocity (rad/s)');
nexttile; plot(c.log_V_MN_ext.Time,c.log_V_MN_ext.Data);hold on;plot(c.log_V_MN_flex.Time,c.log_V_MN_flex.Data,'--');ylabel('MN voltage (mV)');legend('Extensor','Flexor');
nexttile; plot(c.log_A_ext.Time,c.log_A_ext.Data);hold on;plot(c.log_A_flex.Time,c.log_A_flex.Data,'--');ylabel('Activation');xlabel('Time (s)');
for ax=findall(fig,'Type','axes')', grid(ax,'on'); xline(ax,25,':'); xline(ax,30,':'); end
exportgraphics(fig,'settling.pdf','ContentType','vector');exportgraphics(fig,'settling.png','Resolution',150);savefig(fig,'settling.fig');close(fig);
root=fileparts(fileparts(here)); old=load(fullfile(root,'results','sns_demo_results.mat')); twin=load(fullfile(root,'results','knee_twin_ref.mat'));
fprintf('Existing original result: theta %.6g to %.6g deg\n',old.out.log_th.Data(1)*180/pi,old.out.log_th.Data(end)*180/pi);
fprintf('Existing twin: theta %.6g to %.6g deg\n',twin.theta(1)*180/pi,twin.theta(end)*180/pi);
diary off;
end
