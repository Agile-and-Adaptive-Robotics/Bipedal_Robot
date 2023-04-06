%Copy existing 4 figures for isometric torque in the pinned knee joint,
%extensor configuration. Combine into one figure with four tiles.
clear;
clc;
close all;

%First Figure
h1 = openfig('Ext10mm_pin_48cm.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
%Second Figure
h2 = openfig('Ext10mm_pin_46cm.fig','reuse');
ax2 = gca;
%Third Figure
h3 = openfig('Ext10mm_pin_42cm.fig','reuse');
ax3 = gca;
%Fourth Figure
h4 = openfig('Ext10mm_pin_42cm_tendon.fig','reuse');
ax4 = gca;

%Properties to copy
% prop = {'children','FontSize','FontWeight','LineWidth','FontName','XLim','YLim'};

fig5 = figure;
fig5.Position = [197,93,851,669];
h5 = tiledlayout(2,2); %create new figure

s1 = nexttile(1); %create and get handle to the subplot axes
fig1 = get(ax1,'children'); %get handle to all the children in the figure
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
title(s1,'\bf $l_{rest}=48$ cm','Interpreter','latex')
set(s1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 35],'XTick',[-120 -90 -60 -30 0 30],'YLim',[0 14])
lgd1 = legend;
g1 = get(s1,'Children');
set(s1,'Children',[g1(3) g1(2) g1(1)]);

s2 = nexttile(2);
fig2 = get(ax2,'children');
copyobj(fig2,s2);
title(s2,'\bf $l_{rest}=45.7$ cm','Interpreter','latex')
set(s2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 35],'XTick',[-120 -90 -60 -30 0 30],'YLim',[0 14])
lgd2 = legend;
g2 = get(s2,'Children');
set(s2,'Children',[g2(3) g2(2) g2(1)]);

s3 = nexttile(3); %create and get handle to the subplot axes
fig3 = get(ax3,'children');
copyobj(fig3,s3);
title(s3,'\bf $l_{rest}=41.5$ cm, no tendon','Interpreter','latex')
set(s3,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 35],'XTick',[-120 -90 -60 -30 0 30],'YLim',[0 14])
lgd3 = legend('Location','SouthWest');

s4 = nexttile(4);
fig4 = get(ax4,'children');
copyobj(fig4,s4);
title(s4,'\bf $l_{rest}=41.5$ cm, $38$ mm tendon','Interpreter','latex')
set(s4,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 35],'XTick',[-120 -90 -60 -30 0 30],'YLim',[0 14])
lgd4 = legend;

title(h5,'\bf Isometric Torque vs. \theta_{k}, Pinned Knee, Extensor','Interpreter','tex')
ylabel(s1,'\bf Torque, N{\cdot}m','Interpreter','tex')
ylabel(s3,'\bf Torque, N{\cdot}m','Interpreter','tex')
xlabel(s3,'\bf Knee Angle, \circ','Interpreter','tex')
xlabel(s4,'\bf Knee Angle, \circ','Interpreter','tex')
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
    [0.509048178613395 0.87593423019432 0.0366627497062278 0.0573991031390134],...
    'String','\bf (B)',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig5,'textbox',...
    [0.0390129259694477 0.415545590433483 0.036662749706228 0.0573991031390135],...
    'String','\bf (C)',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(fig5,'textbox',...
    [0.50246768507638 0.416741405082213 0.0366627497062278 0.0573991031390135],...
    'String','\bf (D)',...
    'FontSize',12,...
    'FontName','Arial',...
    'FitBoxToText','off',...
    'EdgeColor','none');
