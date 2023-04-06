%Copy existing 4 figures for isometric torque in the pinned knee joint,
%extensor configuration. Combine into one figure with four tiles.
clear;
clc;
close all;

%First Figure
h1 = openfig('Fmax_3D_fit.fig','reuse'); % open figure
f1 = gcf;
% t1 = tiledlayout(h1,1,1);
ax1 = gca; % get handle to axes of figure
cb = colorbar(ax1);
cm = ax1.Colormap;
c1 = findobj(f1, 'type', 'axes');
v1 = c1.View;
fig1 = get(ax1,'children'); %get handle to all the children in the figure
%Second Figure
h2 = openfig('Fmax_620kPa.fig','reuse');
f2 = gcf;
ax2 = gca;
fig2 = get(ax2,'children');
%Third Figure
h3 = openfig('MaxStrainVsRestLength.fig','reuse');
f3 = gcf;
ax3 = gca;
fig3 = get(ax3,'children');


%Properties to copy
% prop = {'children','FontSize','FontWeight','LineWidth','FontName','XLim','YLim'};

fig5 = figure;
fig5.Position = [197,93,851,669];
h5 = tiledlayout(fig5,3,2); %create new figure

s1 = nexttile(1,[2 2]); %create and get handle to the subplot axes
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
title(s1,'\bf $F_{max_{10}}$ vs. $l_{rest}$ and $P$','Interpreter','latex')
set(s1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[0 1],'YLim',[0 800],'ZLim',[0 500]);
lgd1 = legend;

s2 = nexttile(5);
copyobj(fig2,s2);
title(s2,'\bf $F_{max_{10}}$ vs. $l_{rest}$ at P=620 kPa','Interpreter','latex')
set(s2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[0 1],'YLim',[0 500])
lgd2 = legend;


s3 = nexttile(6); %create and get handle to the subplot axes
copyobj(fig3,s3);
title(s3,'\bf \epsilon_{max} vs. l_{rest}','Interpreter','tex')
set(s3,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[0 1],'YLim',[0.1 0.19])
lgd3 = legend('Location','SouthWest');

title(h5,'Max. Force and Strain verus Resting Length')
zlabel(s1,'F_{max_{10}}, N','Interpreter','tex')
xlabel(s1,'l_{rest}, m','Interpreter','tex')
ylabel(s1,'Pressure, kPa','Interpreter','tex')

ylabel(s2,'F_{max_{10}}, N','Interpreter','tex')
xlabel(s2,'l_{rest}, m','Interpreter','tex')

ylabel(s3,'\bf \epsilon_{max}','Interpreter','tex')
xlabel(s3,'\bf l_{rest}, m','Interpreter','tex')

s1.View = v1;
s1.Colormap = cm;

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
