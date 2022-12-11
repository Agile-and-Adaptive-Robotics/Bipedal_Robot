%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
close all;

load KneeExtPin_10mm_all.mat
Theoretical = Torque_46cm';

rest = 0.457;      %resting length, m
kmax = 0.380;               %Length at maximum contraction, m

%% Test 1 done with CALT load cell. Tests 2 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet ExtTest10mm_2 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet ExtTest10mm_1 from Results_table10mm_pinned_FishScale
%% Torque calculated from measurements

Angle1 = [-125	-120.5	-115	-111.5	-108	-101.5	-88	-80	-71	-58	-40	-32	-17.5	-10.5	-4	2	-3.5	-16	-33.5	-48	-58.5	-72	-80	-85	-88.5	-99.5	-111.5	-114.5	-123];
Angle2 = [-102	-100	-81.5	-66	-47	-35	-17	-11	-5];
Angle = [Angle1, Angle2];

Torque1 = [6.168692258	6.072984004	5.91673494	5.983690138	6.220291194	5.580594946	4.981638504	4.559983227	4.059476693	3.15664206	2.56531606	2.127108003	1.785781666	0.774715989	0.309696885	-0.057265164	0.536074229	2.318668963	2.990337108	3.61356103	3.999580201	4.737782449	5.291362695	5.556525636	6.028604398	6.318428803	6.41672146	6.668472107	6.512037525];
Torque2 = [5.969270431	5.781894229	4.978853364	3.667219951	2.623266826	2.087906249	1.579313702	1.311633413	0];
Torque = [Torque1, Torque2];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength1 = [453.5	452	456	450.5	449.5	442.5	439	435	422.5	420	410.5	398.5	400.5	395	384	381	382	394	408.5	414.5	420.5	425	427	432	438	438	441	445	447]/1000;
InflatedLength2 = [430	430	420	415	410	401	390	388	385]/1000;
InflatedLength = [InflatedLength1, InflatedLength2];

ICRtoMuscle1 = [29.5	30	30	29.5	30.5	30	30	32.5	34	35	37	38.5	44.5	50	54	60	55	48	34.5	33	32.5	31.5	30	30	30	30	30	30	30]/1000;
ICRtoMuscle2 = [30	30	30	30	34	35	45	50	55]/1000;
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2];

F1 = zeros(1,size(InflatedLength1, 2));
F2 = zeros(1,size(InflatedLength2, 2));
F = [F1, F2];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand = [TorqueHand1, TorqueHand2];

%load pressure where applicable
test = 2;
runsperseries = 29;

    pres1 = zeros(1,runsperseries);
    
    for j = 1:runsperseries
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test,j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
    end

pres2 = 606*ones(1,size(InflatedLength2, 2));
pres = [pres1 pres2];

for i = 1:size(InflatedLength, 2)  
    F(1,i,1) = festo3(InflatedLength(i), rest, 10, pres(i), kmax);    
    TorqueHand(i) = ICRtoMuscle(i)*F(i);  %Torque will be positive because it is causing extension
end

KMAX = (rest-kmax)/rest;
rel = ((rest-InflatedLength)/rest)/KMAX;
Fn = bpaForce10(rest,rel,pres);

for i = 2:(size(Fn,3)+1)
    F(:,:,i) = Fn(:,:,i-1);
end

for i = 1:size(F,3)
    TorqueHand(:,:,i) = ICRtoMuscle.*F(1,:,i);  %Torque will be positive because it is causing extension   
    TorqueHand1 = TorqueHand(1,1:size(TorqueHand1,2),i);
    TorqueHand2 = TorqueHand(1,((size(TorqueHand1,2)+1)):size(TorqueHand,2),i); 
end


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
Ma = Vas_Pam_46cm.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
ax11 = subplot(2,1,1);
hold on
title('\bf Expected vs measured $r_{\hat{k}}$', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
pp21 = plot(phiD,G,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected $r_{\hat{k}}$');
ss2 = scatter(Angle, ICRtoMuscle,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured $r_{\hat{k}}$');
set(ax11,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdMa11 = legend('Interpreter','latex');
lgdMa11.FontSize = 12;
hold off

ax12 = subplot(2,1,2);
hold on
title('\bf Expected vs measured $r_{\hat{k}}$', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
pp22 = plot(phiD,G,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected $r_{\hat{k}}$');
ss2_1 = scatter(Angle1, ICRtoMuscle1,sz,'filled','MarkerFaceColor',c1,'DisplayName','\bf Measured $r_{\hat{k}}$, BB{\&}JM LC');
ss2_2 = scatter(Angle2, ICRtoMuscle2,sz,'filled','MarkerFaceColor',c2,'DisplayName','\bf Measured $r_{\hat{k}}$, BB FS');
set(ax12,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial')
lgdMa12 = legend('Interpreter','latex');
lgdMa12.FontSize = 12;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = (rest-(Vas_Pam_46cm.MuscleLength-tendon-2*fitting))/rest;
relstrain = (strain)./KMAX;
realRel = (rest-InflatedLength)/rest/KMAX;
realRel1 = (rest-InflatedLength1)/rest/KMAX;
realRel2 = (rest-InflatedLength2)/rest/KMAX;

figure
ax21 = subplot(2,1,1);
hold on
title('Expected vs measured \epsilon^*','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
plot(phiD,relstrain,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected \epsilon^*')
scatter(Angle,realRel,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured \epsilon^*')
set(ax21,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdEp1 = legend('Interpreter','tex');
lgdEp1.FontSize = 12;
hold off

ax22 = subplot(2,1,2);
hold on
title('Expected vs measured \epsilon^*','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
plot(phiD,relstrain,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected \epsilon^*')
scatter(Angle1,realRel1,sz,'filled','MarkerFaceColor',c1,'DisplayName','\bf Measured \epsilon^*, BB+JM LC')
scatter(Angle2,realRel2,sz,'filled','MarkerFaceColor',c2,'DisplayName','\bf Measured \epsilon^*, BB FS')
set(ax22,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdEp2 = legend('Interpreter','tex');
lgdEp2.FontSize = 12;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Vas_Pam_46cm.MuscleLength-2*fitting-tendon;

figure
hold on
title('\bf Expected vs measured {l_{m}}','FontWeight','bold','Interpreter','tex')
xlabel('\bf Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('\bf {l_{m}}, m','FontWeight','bold','Interpreter','tex')
plot(phiD,MuscleLength,'DisplayName','\bf Expected {l_{m}}')
scatter(Angle,InflatedLength,'DisplayName','\bf Measured {l_{m}}')
axLm = gca;
set(axLm,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdLm = legend('Interpreter','tex');
lgdLm.FontSize = 12;
hold off

%% Plotting Z axis torque values

figure
hold on
title('Iso. Torque vs {\theta_{k}}, Pinned, Extensor, l_{rest} = 45.7cm','Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
gca1 = gca;
gcf1 = gcf;

PL1 = plot(phiD, Theoretical,'Color',c4,'Linewidth',2,'DisplayName','Theoretical Torque');
scM = scatter(Angle,Torque,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Measured Torque');
scH = scatter(Angle,TorqueHand(:,:,4),sz2,'filled','MarkerFaceColor',c1,'DisplayName','Back calculated Torque');

set(gca1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial')
lgd1 = legend;
lgd1.FontSize = 12;
hold off

%% Plotting Z axis torque values, color coded to mark separate tests

figure
hold on
title('Iso. Torque vs {\theta_{k}}, Pinned, Extensor, l_{rest} = 45.7cm','Interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','Interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','Interpreter','tex')
gca2 = gca;
gcf2 = gcf;

PL2 = plot(phiD, Theoretical,'Color',c7,'Linewidth',2,'DisplayName','Theoretical');
scM1 = scatter(Angle1,Torque1,sz,'d','filled','MarkerFaceColor',c6,'DisplayName','\bf BB{\&}JM LC, measured');
scH1 = scatter(Angle1,TorqueHand1,sz2,'filled','MarkerFaceColor',c4,'DisplayName','\bf BB{\&}JM LC, back calc');
scM2 = scatter(Angle2,Torque2,sz,'d','filled','MarkerFaceColor',c2,'DisplayName','\bf BB FS, measured');
scH2 = scatter(Angle2,TorqueHand2,sz2,'filled','MarkerFaceColor',c1,'DisplayName','\bf BB FS, back calc');

set(gca2,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial')
lgd2 = legend('Interpreter','latex');
lgd2.FontSize = 12;
hold off

%% Mean and RMSE
Tqz = cell(1,1);
Tqz{1} = Vas_Pam_46cm.Torque(:,3,4);        %Calculated Torque, new simplified exponential equation w/o optimized fitting length
%Tqz{2} = Vas_Pam_46cm_adj.Torque(:,3,4);   %Calculated Torque, adjusted with optimized fitting length
%Tqz{3} = TorqueHand(:,:,4);                %Placeholder in case we want to compare SSE/RMSE of back calculated torque to measured torque

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