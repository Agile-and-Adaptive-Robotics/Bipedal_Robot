%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
close all;

% f = fullfile('github/bipedal_robot/code/matlab');
% qt = addpath(genpath(f));

load KneeFlx_10mm_42cm.mat
Theoretical = TorqueR(:,3)';
%rest = 0.415, tendon = 0.012

%% Optimized Theoretical Calculation
Theo_adj = TorqueR_adj(:,3)';

% rest = 0.415; %set resting length back to measured value
% tendon = 0.012; %set tendon back to measured(?) value
% kmax = 0.350; %set kmax back to measured value
%% Test 1 done with CALT load cell
%Test 1 == sheet FlxTest10mm from Results_table_10mm

%% Torque calculated from measurements
Angle = [-2	-11	-16	-26.5	-40	-49.5	-54.5	-61	-63.5	-81];

Torque = [-8.067539602	-7.933051945	-7.313409814	-5.440386864	-3.860259596	-2.60401704	-2.31467059	-1.34002924	-0.8951295	-0.355778987];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength = [397	390	387	380	375	371	367	365	365	357]/1000;

ICRtoMuscle = [38	42	42	40	35	39	35	35	33	29]/1000;

F = zeros(size(InflatedLength));

%load pressure where applicable
runsperseries = 10;

pres = zeros(1,runsperseries);

     for j = 1:runsperseries
                file_name = sprintf('FlxTest%0.0f.mat', j);
                load(file_name,'Stats')
                pres(1,j) = Stats{'Mean',2};
     end

KMAX = (rest-kmax)/rest;
rel = ((rest-InflatedLength)/rest)/KMAX;     
     
for i = 1:length(InflatedLength)  
    F(i) = bpaForce10(rest, rel(i), pres(i));    
end

TorqueHand = -ICRtoMuscle.*F;  %Torque will be negative because it is causing flexion

%% Plot expected versus measured moment arm
Ma = Bifemsh_Pam_adj.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
hold on
pp2 = plot(phiD,G,'DisplayName','\bf Expected $r_{\hat{k}}$');
ss2 = scatter(Angle, ICRtoMuscle,'DisplayName','\bf Measured $r_{\hat{k}}$');
title('\bf Expected vs measured $r_{\hat{k}}$', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 12;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','latex');
lgdMa.FontSize = 12;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = (rest-(Bifemsh_Pam_adj.MuscleLength-tendon-2*fitting))/rest;
relstrain = (strain)./KMAX;
realRel = (rest-InflatedLength)/rest/KMAX;

figure
hold on
plot(phiD,relstrain,'DisplayName','\bf Expected \epsilon^*')
scatter(Angle,realRel,'DisplayName','\bf Measured \epsilon^*')
title('Expected vs measured \epsilon^*','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
ax2 = gca;
ax2.FontSize = 12;
ax2.FontWeight = 'bold';
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 12;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','tex');
lgdMa.FontSize = 12;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Bifemsh_Pam_adj.MuscleLength-2*fitting-tendon;

figure
hold on
plot(phiD,MuscleLength,'DisplayName','\bf Expected $l_{m}$')
scatter(Angle,InflatedLength,'DisplayName','\bf Measured $l_{m}$')
title('\bf Expected vs measured $l_{m}$','Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf $l_{m}$, m','Interpreter','latex')
ax3 = gca;
ax3.FontSize = 12;
ax3.FontWeight = 'bold';
ax3.FontName = 'Arial';
ax3.YAxis.LineWidth = 2; ax3.YAxis.FontSize = 12;
ax3.XAxis.LineWidth = 2; ax3.XAxis.FontSize = 12;
lgdMa = legend('Interpreter','latex');
lgdMa.FontSize = 12;
hold off

%% Plotting with polynomial solver
%Matlab hex color values:
c1 = '#FFD700'; %gold
c2 = '#FFB14E'; %orange
c3 = '#FA8775'; %light orange
c4 = '#EA5F94'; %pink
c5 = '#CD34B5'; %magenta
c6 = '#9D02D7'; %magenta 2
c7 = '#0000FF'; %indigo
c8 = '#000000'; %black
sz = 60;        %size of data points
c = {c1; c2; c3; c4; c5; c6; c7; c8};

figure('units','normalized')
hold on
title('\bf Torque for Biomimetic Knee, $10\,$mm Flexor, $l_{rest}=41.5\,$cm','interpreter','latex')
xlabel('Knee angle, \circ','FontWeight','bold','interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','interpreter','tex')
gca1 = gca;
gcf1 = gcf;
set(gca,'FontSize', 12, 'FontWeight', 'bold');
PL1 = plot(phiD, Theoretical,'Color',c1,'Linewidth',2,'DisplayName','Theoretical, original');
chr = 'Theoretical, optimized';
PL3 = plot(phiD, Theo_adj,'Color',c3,'Linewidth',2,'DisplayName',chr);
sc5 = scatter(Angle,TorqueHand,sz,'filled','MarkerFaceColor',c5,'DisplayName','Back calculated torque');
sc7 = scatter(Angle,Torque,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Measured Torque');
ax4 = gca;
ax4.FontSize = 12;
ax4.FontWeight = 'bold';
ax4.FontName = 'Arial';
ax4.YAxis.LineWidth = 2; ax4.YAxis.FontSize = 10;
ax4.XAxis.LineWidth = 2; ax4.XAxis.FontSize = 10;
%lgd = legend('interpreter','latex');
lgd = legend;
hold off

%% Mean and RMSE
% Tqz = cell(2,1);
% Tqz{1} = Bifemsh_Pam.Torque(:,3,4);         %Calculated Torque, new simplified exponential equation w/o optimized fitting length
% Tqz{2} = Bifemsh_Pam_adj.Torque(:,3,4);     %Calculated Torque, adjusted with optimized fitting length
% %Tqz{3} = TorqueHand(:,:,4);                %Placeholder in case we want to compare SSE/RMSE of back calculated torque to measured torque
% 
% %fit options
% mod_Pam = fittype('cubicinterp');
% Options = fitoptions(mod_Pam);
% Options.Normal = 'on';
% 
% %prepare cells
% mdl_Pam = cell(size(Tqz,1));
% val = cell(length(mdl_Pam),1);
%           
% %Get values at each angle there is measurement data for
% for j = 1:length(Tqz)
%      Options.Exclude = isnan(Tqz{j});
%      mdl_Pam{j} = fit(phiD',Tqz{j},mod_Pam,Options);
%      val{j} = feval(mdl_Pam{j},Angle');
% end
% 
% y = Torque';        
% ynew = val;
%         
% yresid = cell(length(ynew),1);
% SSresid = cell(length(ynew),1);
% fu = cell(length(ynew),1);
%         
% for i = 1:length(ynew)
%     yresid{i} = y-ynew{i};              %residual error
%     SSresid{i} = sum(yresid{i}.^2,'omitnan'); %Sum of squares of the residual
%     fu{i} = sqrt(SSresid{i}/length(yresid{i}));        % RMSE for function 1
% end
% 
% fprintf('Original torque calculation returns SSE of %5d with an RMSE of %5d\n',SSresid{1},fu{1})
% fprintf('Optimized torque calculation returns SSE of %5d with an RMSE of %5d\n',SSresid{2},fu{2})