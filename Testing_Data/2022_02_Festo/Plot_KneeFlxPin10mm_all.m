%Copy existing 4 figures for the pinned knee joint,
%flexor configuration. Combine into one figure with four tiles.
clear;
clc;
close all;

%First Figure
h1 = openfig('Flx10mm_pin_48cm.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
%Second Figure
h2 = openfig('Flx10mm_pin_48cm_relStrain.fig','reuse');
ax2 = gca;
%Third Figure
h3 = openfig('Flx10mm_pin_48cm_optimized.fig','reuse');
ax3 = gca;
ylabel('Relative strain \epsilon^*');
%Fourth Figure
h4 = openfig('Flx10mm_pin_46cm_optimized.fig','reuse');
ax4 = gca;

xLim = [-120 15];
xTik = [-120 -90 -60 -30 0];

fig5 = figure;
fig5.Position = [197,93,851,669];
h5 = tiledlayout(2,2); %create new figure

s1 = nexttile(1); %create and get handle to the subplot axes
fig1 = get(ax1,'children'); %get handle to all the children in the figure
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
% title(s1,'\bf l_{rest}=48.5 cm')
set(s1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',xLim,'XTick',xTik)
set(s1,'XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05])
lgd1 = legend('Location','SouthWest','FontSize',8);

s2 = nexttile(2);
fig2 = get(ax2,'children');
copyobj(fig2,s2);
% title(s2,'\bf l_{rest}=48.5 cm')
set(s2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',xLim,'XTick',xTik)
set(s2,'XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05])
lgd2 = legend('Location','SouthWest','FontSize',8);

s3 = nexttile(3); %create and get handle to the subplot axes
fig3 = get(ax3,'children');
copyobj(fig3,s3);
% title(s3,'\bf l_{rest}=48.5 cm, optimized')
set(s3,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-120 15],'XTick',xTik,'YLim',[-30 0])
set(s3,'XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05])
lgd3 = legend('Location','SouthWest','FontSize',8);

s4 = nexttile(4);
fig4 = get(ax4,'children');
copyobj(fig4,s4);
% title(s4,'l_{rest}=45.7 cm, validated')
set(s4,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XTick',xTik,'YLim',[-30 0])
set(s4,'XMinorTick','on','YMinorTick','on','TickLength',[0.025, 0.05])
lgd4 = legend('Location','SouthWest','FontSize',8);

title(h5,'\bf Optimization of fitting length')
ylabel(s1,'\bf Torque, N\cdotm','Interpreter','tex')
ylabel(s3,'\bf Torque, N\cdotm','Interpreter','tex')
xlabel(s3,'\bf Knee angle, \circ')
xlabel(s4,'\bf Knee angle, \circ')
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
