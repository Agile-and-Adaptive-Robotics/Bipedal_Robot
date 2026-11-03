% Restore the BFLH two-coordinate calculation specified in the archived2019
% sources. Source code/input MATs are read-only. This does not redesign or fit.
outDir=fileparts(mfilename('fullpath'));repoRoot=outDir;
while true
 [par,nm]=fileparts(repoRoot);if strcmpi(nm,'Bipedal_Robot'),break;end
 assert(~strcmp(par,repoRoot));repoRoot=par;
end
legacy=fullfile(repoRoot,'tmp','old_code_copy_20261022');
addpath(legacy,'-begin');oldFolder=pwd;cleanup=onCleanup(@()cd(oldFolder));cd(legacy);
set(groot,'DefaultFigureVisible','off');warning('off','all');
H=load(fullfile(legacy,'AttachPoints.mat'));
R=load(fullfile(legacy,'AttachPoints_Robot.mat'));
M=load(fullfile(legacy,'MuscProp.mat'));
P=load(fullfile(legacy,'PAM_Prop.mat'));
K=load(fullfile(legacy,'MovePoints.mat'));
N=100;knee=linspace(-120,10,N+1);knee=knee(1:N);
hipHuman=linspace(-20,140,N+1);hipHuman=hipHuman(1:N);
hipRobot=linspace(-15,110,N+1);hipRobot=hipRobot(1:N);
human=zeros(N,N);robot=zeros(N,N);robotLength=zeros(N,N);robotArm=zeros(N,N);p_h=[-.0707;-.0661;.0835];
for j=1:N
 ThH=[rz(deg2rad(hipHuman(j))),p_h;0 0 0 1];
 ThR=[rz(deg2rad(hipRobot(j))),p_h;0 0 0 1];
 for i=1:N
  q=deg2rad(knee(i));
  Tk=[rz(q),[K.fcn1(q);K.fcn2(q);0];0 0 0 1];
  lh=norm(vt(ThH*Tk,H.bifemlh_wr1)-H.bifemlh_o)+norm(H.bifemlh_i-H.bifemlh_wr1);
  dirh=vt((ThH*Tk)\eye(4),H.bifemlh_o)-H.bifemlh_wr1;dirh=dirh/norm(dirh);
  mah=cross(H.bifemlh_wr1,dirh);
  fh=forz(lh,M.bifemlh_mif,M.bifemlh_ofl,M.bifemlh_tsl,M.bifemlh_pa);
  human(i,j)=fh*mah(3);
  lr=norm(vt(ThR*Tk,R.twr2)-R.p5)+norm(R.twr3-R.twr2)+norm(R.twr3-R.t5);
  dirr=vt((ThR*Tk)\eye(4),R.p5)-R.twr2;dirr=dirr/norm(dirr);
  mar=cross(R.twr2,dirr);
  robotLength(i,j)=lr;robotArm(i,j)=mar(3);
 end
 if mod(j,20)==0,fprintf('BFLH %d/%d columns\n',j,N);end
end
% The archived PAM_Prop file supplies MIF; the original script's Size() call
% derives diameter, active rest length, and longest route from these lengths.
evalc('[diameter,restLength,longestRoute]=Size(robotLength,P.bifemlh_mif);');
fprintf('Historical robot sizing: diameter %g mm, rest %g m, longest route %g m\n',diameter,restLength,longestRoute);
for j=1:N
 for i=1:N
  robot(i,j)=festo(robotLength(i,j),restLength,diameter,longestRoute)*robotArm(i,j);
 end
end
assert(isreal(human)&&isreal(robot));
fprintf('Historical BFLH: human [%g,%g], robot [%g,%g] Nm\n',min(human(:)),max(human(:)),min(robot(:)),max(robot(:)));
save(fullfile(outDir,'historical_bflh_surfaces.mat'),'human','robot','knee','hipHuman','hipRobot');
zlo=10*floor(min([human(:);robot(:)])/10);cl=[zlo 0];
fig=figure('Color','w','Position',[20 20 1680 760]);t=tiledlayout(fig,1,2,'TileSpacing','loose','Padding','loose');
for p=1:2
 ax=nexttile(t);if p==1,x=hipHuman;z=human;else,x=hipRobot;z=robot;end
 surf(ax,x,knee,z,'EdgeColor','none');
 set(ax,'FontName','Arial','FontSize',17,'XLim',[-20 140],'YLim',[-120 10], ...
  'ZLim',cl,'CLim',cl,'XTick',[0 50 100],'YTick',[-120 -60 0], ...
  'ZTick',zlo:20:0,'LineWidth',1);
 colormap(ax,parula(256));view(ax,[142.5 30]);grid(ax,'on');
 xlabel(ax,'Hip angle (°)');ylabel(ax,'Knee angle (°)');zlabel(ax,'Knee torque (N·m)');
 cb=colorbar(ax);cb.Limits=cl;cb.Ticks=zlo:20:0;
end
exportgraphics(fig,fullfile(outDir,'historical_bflh_common_scales.png'),'Resolution',180);
exportgraphics(fig,fullfile(outDir,'historical_bflh_common_scales.pdf'),'ContentType','vector');
savefig(fig,fullfile(outDir,'historical_bflh_common_scales.fig'));
fprintf('Identical z/color scales: %s\n',mat2str(cl));close(fig);
function R=rz(q)
 R=[cos(q),-sin(q),0;sin(q),cos(q),0;0 0 1];
end
function p=vt(T,v)
 h=T*[v(:);1];p=h(1:3);
end
