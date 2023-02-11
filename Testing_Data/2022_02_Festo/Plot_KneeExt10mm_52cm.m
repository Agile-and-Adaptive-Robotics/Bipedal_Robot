%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
close all;

load KneeExt_10mm_52cm.mat
rest = 0.520;      %resting length, m
kmax = 0.440;               %Length at maximum contraction, m
tendon = 0.030;             %Tendon, measured
% Vas_Pam = MonoPamDataExplicit(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);
Torque_52cm = TorqueR(:,3);
Theoretical = Torque_52cm';

%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm from Results_table_10mm
%% Torque calculated from measurements

Angle = [-120	-108.5	-93	-70	-53.5	-49	-21	-9	-2];

Torque = TorqueZ;

%% Calculate Torque by finding force from muscle contraction and distance

InflatedLength = [518	520	510	495	485	475	460	450	450]/1000;

ICRtoMuscle = [38	35	40	40	45	45	40	40	41]/1000;



%load pressure where applicable
test = 1;
runsperseries = 9;

    pres = zeros(1,runsperseries);
    
    for j = 1:runsperseries
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres(1,j) = Stats{'Mean',2};
    end

KMAX = (rest-kmax)/rest;
rel = ((rest-InflatedLength)/rest)/KMAX;
Fn = bpaForce10(rest,rel,pres);

F = Fn(3);
TorqueHand = ICRtoMuscle.*F;  %Torque will be positive because it is causing extension   

%% Prepare for plotting
% Create accessible color scheme
c1 = '#FFD700'; %gold
c2 = '#FFB14E'; %orange
c3 = '#FA8775'; %light orange
c4 = '#EA5F94'; %pink
c5 = '#CD34B5'; %magenta
c6 = '#9D02D7'; %magenta 2
c7 = '#0000FF'; %indigo
c8 = '#000000'; %black
sz = 60;        %size of data points
sz2 = sz*0.666; %size of second data points
c = {c1; c2; c3; c4; c5; c6; c7; c8};

%% Plot expected versus measured moment arm
Ma = Vas_Pam.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
ax11 = gca;
hold on
title('\bf Expected vs measured $r_{\hat{k}}$', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
pp21 = plot(phiD,G,'Color',c4,'LineWidth',2,'DisplayName','\bf Expected $r_{\hat{k}}$');
ss2 = scatter(Angle, ICRtoMuscle,sz,'filled','MarkerFaceColor',c7,'DisplayName','\bf Measured $r_{\hat{k}}$');
set(ax11,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdMa11 = legend('Interpreter','latex');
lgdMa11.FontSize = 12;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = (rest-(Vas_Pam.MuscleLength-tendon-2*fitting))/rest;
relstrain = (strain)./KMAX;
realRel = (rest-InflatedLength)/rest/KMAX;


figure
ax21 = gca;
hold on
title('Expected vs measured \epsilon^*','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
plot(phiD,relstrain,'Color',c4,'LineWidth',2,'DisplayName','\bf Expected \epsilon^*')
scatter(Angle,realRel,sz,'filled','MarkerFaceColor',c7,'DisplayName','\bf Measured \epsilon^*')
set(ax21,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdEp1 = legend('Interpreter','tex');
lgdEp1.FontSize = 12;
hold off


%% Plot measured versus expected BPA length
MuscleLength = Vas_Pam.MuscleLength-2*fitting-tendon;

figure
hold on
title('\bf Expected vs measured {l_{m}}','FontWeight','bold','Interpreter','tex')
xlabel('\bf Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('\bf {l_{m}}, m','FontWeight','bold','Interpreter','tex')
plot(phiD,MuscleLength,'Color',c4,'DisplayName','\bf Expected {l_{m}}')
scatter(Angle,InflatedLength,'MarkerFaceColor',c7,'DisplayName','\bf Measured {l_{m}}')
axLm = gca;
set(axLm,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdLm = legend('Interpreter','tex');
lgdLm.FontSize = 12;
hold off

%% Plotting Z axis torque values

figure
hold on
title('Iso. Torque vs {\theta_{k}}, Pinned, Extensor, l_{rest} = 52.0 cm','Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
gca1 = gca;
gcf1 = gcf;

PL1 = plot(phiD, Theoretical,'Color',c4,'Linewidth',2,'DisplayName','Theoretical');
scM = scatter(Angle,Torque,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Measured');
scH = scatter(Angle,TorqueHand,sz2,'filled','MarkerFaceColor',c1,'DisplayName','Back calculated');

set(gca1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 40])
lgd1 = legend;
lgd1.FontSize = 12;
hold off

%% Mean and RMSE
Tqz = cell(1,1);
Tqz{1} = Vas_Pam.Torque(:,3);        %Calculated Torque, new simplified exponential equation w/o optimized fitting length
%Tqz{2} = Vas_Pam_adj.Torque(:,3);   %Calculated Torque, adjusted with optimized fitting length
%Tqz{3} = TorqueHand;                %Placeholder in case we want to compare SSE/RMSE of back calculated torque to measured torque

%fit options
mod_Pam = fittype('cubicinterp');
Options = fitoptions(mod_Pam);
Options.Normal = 'on';

%prepare cells
mdl_Pam = cell(size(Tqz,1));
val = cell(length(mdl_Pam),1);
          
%Get values at each angle there is measurement data for
for j = 1:length(Tqz)
     Options.Exclude = isnan(Tqz{j});
     mdl_Pam{j} = fit(phiD',Tqz{j},mod_Pam,Options);
     val{j} = feval(mdl_Pam{j},Angle');
end

y = Torque';        
ynew = val;
        
yresid = cell(length(ynew),1);
SSresid = cell(length(ynew),1);
fu = cell(length(ynew),1);
        
for i = 1:length(ynew)
    yresid{i} = y-ynew{i};              %residual error
    SSresid{i} = sum(yresid{i}.^2,'omitnan'); %Sum of squares of the residual
    fu{i} = sqrt(SSresid{i}/length(yresid{i}));        % RMSE for function 1
end

fprintf('Original torque calculation returns SSE of %5d with an RMSE of %5d\n',SSresid{1},fu{1})