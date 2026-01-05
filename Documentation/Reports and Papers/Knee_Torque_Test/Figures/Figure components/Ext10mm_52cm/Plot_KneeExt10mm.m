%Copy existing 4 figures for isometric torque in the pinned knee joint,
%extensor configuration. Combine into one figure with four tiles.
clear;
clc;
close all;

%First Figure
h1 = openfig('Ext10mm_52cm.eps','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
%Second Figure
h2 = openfig('Ext10mm_52cm_MomentArm.eps','reuse');
ax2 = gca;
%Third Figure
h3 = openfig('Ext10mm_52cm_strain.fig','reuse');
ax3 = gca;


%Properties to copy
% prop = {'children','FontSize','FontWeight','LineWidth','FontName','XLim','YLim'};

fig5 = figure;
fig5.Position = [197,93,851,669];
h5 = tiledlayout(2,2); %create new figure

s1 = nexttile(1); %create and get handle to the subplot axes
fig1 = get(ax1,'children'); %get handle to all the children in the figure
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
% title(s1,'\bf l_{rest}=48.0 cm','Interpreter','tex')
% set(s1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 35],'XTick',[-120 -90 -60 -30 0 30])
% set(s1,'YLim',[0 15], 'XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
lgd1 = legend('FontSize',8);
g1 = get(s1,'Children');
set(s1,'Children',[g1(3) g1(2) g1(1)]);

s2 = nexttile(2);
fig2 = get(ax2,'children');
copyobj(fig2,s2);
% title(s2,'\bf l_{rest}=45.7 cm','Interpreter','tex')
% set(s2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 35],'XTick',[-120 -90 -60 -30 0 30],'YLim',[0 14])
% set(s2,'YLim',[0 15], 'XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
lgd2 = legend('FontSize',8);
g2 = get(s2,'Children');
set(s2,'Children',[g2(3) g2(2) g2(1)]);

s3 = nexttile(3); %create and get handle to the subplot axes
fig3 = get(ax3,'children');
copyobj(fig3,s3);
% title(s3,'\bf l_{rest}=41.5 cm, no tendon','Interpreter','tex')
% set(s3,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 35],'XTick',[-120 -90 -60 -30 0 30],'YLim',[0 14])
% set(s3,'YLim',[0 15], 'XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05]);
lgd3 = legend('Location','SouthWest','FontSize',8);


title(h5,'\bf Biomimetic Knee, \phi10 mm Extensor','Interpreter','tex')
ylabel(s1,'\bf Torque, N{\cdot}m','Interpreter','tex')
ylabel(s2,'\bf Moment arm, N{\cdot}m','Interpreter','tex')
ylabel(s2,'\bf Relative strain, \epsilon^*','Interpreter','tex')
xlabel(s1,'\bf Knee Angle, \circ','Interpreter','tex')
xlabel(s2,'\bf Knee Angle, \circ','Interpreter','tex')
xlabel(s3,'\bf Knee Angle, \circ','Interpreter','tex')
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
