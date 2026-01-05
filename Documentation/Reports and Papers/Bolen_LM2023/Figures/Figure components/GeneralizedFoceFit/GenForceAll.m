%Copy existing 4 figures for isometric torque in the pinned knee joint,
%extensor configuration. Combine into one figure with four tiles.
clear;
clc;
close all;

%First Figure
h1 = openfig('GeneralizedForceFit10.fig','reuse'); % open figure
f1 = gcf;
% t1 = tiledlayout(h1,1,1);
ax1 = gca; % get handle to axes of figure
cb1 = colorbar(ax1);
cm1 = ax1.Colormap;
c1 = findobj(f1, 'type', 'axes');
v1 = c1.View;

%Second Figure
h2 = openfig('GeneralizedForceFit20.fig','reuse');
f2 = gcf;
ax2 = gca;
cb2 = colorbar(ax2);
cm2 = ax2.Colormap;
c2 = findobj(f2, 'type', 'axes');
v2 = c2.View;

%Third Figure
h3 = openfig('GeneralizedForceFit40.fig','reuse');
f3 = gcf;
ax3 = gca;
cb3 = colorbar(ax3);
cm3 = ax3.Colormap;
c3 = findobj(f3, 'type', 'axes');
v3 = c3.View;



%Properties to copy
% prop = {'children','FontSize','FontWeight','LineWidth','FontName','XLim','YLim'};

fig5 = figure;
fig5.Position = [197,93,851,669];
h5 = tiledlayout(fig5,3,1); %create new figure

s1 = nexttile(1); %create and get handle to the subplot axes
fig1 = get(ax1,'children'); %get handle to all the children in the figure
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
title(s1,'\bf Normalized force for {\phi}10 mm BPAs ','Interpreter','tex')
set(s1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[0 1],'YLim',[0 1],'ZLim',[0 1.1]);
% lgd1 = legend;
s1.View = v1;
s1.Colormap = cm1;

s2 = nexttile(2);
fig2 = get(ax2,'children');
copyobj(fig2,s2);
title(s2,'\bf Normalized force for {\phi}20 mm BPAs ','Interpreter','tex')
set(s2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[0 1],'YLim',[0 1],'ZLim',[0 1.1])

% lgd2 = legend;


s3 = nexttile(3); %create and get handle to the subplot axes
fig3 = get(ax3,'children');
copyobj(fig3,s3);
title(s3,'\bf Normalized force for {\phi}40 mm BPAs ','Interpreter','tex')
set(s3,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[0 1],'YLim',[0 1],'ZLim',[0 1.1])
s2.View = v2;
s2.Colormap = cm2;
% lgd3 = legend('Location','southeast');

zlabel(s1,'F^{\ast}','Interpreter','tex')
xlabel(s1,'\epsilon^{\ast}','Interpreter','tex')
ylabel(s1,'P^{\ast}','Interpreter','tex');

zlabel(s2,'F^{\ast}','Interpreter','tex')
xlabel(s2,'\epsilon^{\ast}','Interpreter','tex')
ylabel(s2,'P^{\ast}','Interpreter','tex');

zlabel(s3,'F^{\ast}','Interpreter','tex')
xlabel(s3,'\epsilon^{\ast}','Interpreter','tex')
ylabel(s3,'P^{\ast}','Interpreter','tex');





s3.View = v3;
s3.Colormap = cm3;

% Create textbox
annotation(fig5,'textbox',...
    [0.036662749706228 0.871748878923767 0.036662749706228 0.0573991031390134],...
    'String','\bf (A)',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');


% Create textbox
annotation(fig5,'textbox',...
    [0.0390129259694477 0.415545590433483 0.036662749706228 0.0573991031390135],...
    'String','\bf (B)',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig5,'textbox',...
    [0.50246768507638 0.416741405082213 0.0366627497062278 0.0573991031390135],...
    'String','\bf (C)',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');
