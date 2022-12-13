%% Pinned knee, Extensor
%Run and save data from testing results
clear;
clc;
%  close all;

load KneeExtPin_10mm_all.mat
rest = 0.480;      %resting length, m
kmax = 0.398;      %Length at maximum contraction, m
Vas_Pam_48cm = MonoPamDataExplicit_compare(Name, Location, CrossPoint, Dia, T, rest, kmax, tendon, fitting, pres);
Torque_48cm = Vas_Pam_48cm.Torque(:,3,4);
Theoretical = Torque_48cm';



%% Tests 1 and 4 done with CALT load cell. Tests 2 and 3 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 


%% Torque calculated from measurements

Angle1 = [-121.5	-108.5	-101.5	-94	-83	-75	-61.5	-65	-32	-18	-32.5	-47.5	-57.5	-69.5	-80	-90.5	-96.5	-104	-113.5	-123];
Angle2 = [-110	-82	-71.5	-55	-42	-43	-65.5	-74.5	-72	-82	-88.5	-98	-113.5];
Angle3 = [-110	-105	-103.5	-91	-80	-72.5	-62	-48	-40.5	-31	-30	-31.5	-43	-50	-62.5	-69	-78	-88	-93	-100	-106	-107];
Angle4 = [-123	-119	-106	-92];
Angle = [Angle1, Angle2, Angle3, Angle4];

Torque1 = [4.62685036	4.503204046	4.29202985	4.002027833	3.417014674	2.826201474	1.538232508	1.214084497	0.439203694	-0.080856656	0.539217698	1.260349872	2.254964352	3.207688111	3.851548806	4.243514283	4.672306017	4.873564883	5.146026298	6.087161098];
Torque2 = [5.783034531	4.71210221	3.667943198	1.847358253	1.044159012	1.579625173	2.356051105	3.159250345	3.50730335	4.203409358	4.203409358	4.97983529	4.97983529];
Torque3 = [5.783034531	4.71210221	4.471142438	3.935676278	3.667943198	3.159250345	2.623784185	1.847358253	1.419266981	0.53546616	0	0.53546616	1.044159012	2.356051105	2.891517265	3.400210118	3.935676278	4.203409358	4.71210221	4.97983529	6.291727383	6.559460463];
Torque4 = [7.013441569	4.72935096	4.346134571	3.454623635];
Torque = [Torque1, Torque2, Torque3, Torque4];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for tests 1 & 2 were done incorrectly and should be disregarded

InflatedLength1 = [455.5	461	451	451	440	438.5	435	415.5	413	410	409.5	418	426.5	438	442.5	441.5	450	458.5	455.5	463]/1000;
InflatedLength2 = [448	445.5	429.5	413.5	401	410	417.5	422.5	430.5	429.5	441.5	439	443.5]/1000;
InflatedLength3 = [445	445	433	426	423	421	420	415	405	400	398	400	405	407	416	420	430	433	441	440	450	453]/1000;
InflatedLength4 = [459	456	447	438]/1000;
InflatedLengthX = [InflatedLength3, InflatedLength4];
InflatedLength = [InflatedLength1, InflatedLength2, InflatedLength3, InflatedLength4];

ICRtoMuscle1 = [30.5	41.5	39	35.5	37	40.5	37	46.5	38.5	45.5	36.5	38.5	35	38.5	36	42	39	36	38	37]/1000;
ICRtoMuscle2 = [46	54.5	53	52	57	49	53	53	54	51	52	53	55.5]/1000;
ICRtoMuscle3 = [32	30	30	30	30	30	30	30	35	35	40	35	34	33	30	30	30	30	30	30	30	30]/1000;
ICRtoMuscle4 = [30	30	30	30]/1000;
ICRtoMuscleX = [ICRtoMuscle3, ICRtoMuscle4];
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2, ICRtoMuscle3, ICRtoMuscle4];

F1 = zeros(1,size(InflatedLength1, 2));
F2 = zeros(1,size(InflatedLength2, 2));
F3 = zeros(1,size(InflatedLength3, 2));
F4 = zeros(1,size(InflatedLength4, 2));
% F = [F1, F3, F4];
F = [F1, F2, F3, F4];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand3 = zeros(1,size(InflatedLength3, 2));
TorqueHand4 = zeros(1,size(InflatedLength4, 2));
% TorqueHand = [TorqueHand1, TorqueHand3, TorqueHand4];
TorqueHand = [TorqueHand1, TorqueHand2, TorqueHand3, TorqueHand4];

%load pressure where applicable
test = [1 3];
runsperseries = [20 4];

    pres1 = zeros(1,size(runsperseries(1),2));
    pres4 = zeros(1,size(runsperseries(2),2));
    
    for j = 1:runsperseries(1)
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(1),j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
    end
    for j = 1:runsperseries(2)
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(2),j);
                load(file_name,'Stats')
                pres4(1,j) = Stats{'Mean',2};
    end
pres2 = 612*ones(1,size(InflatedLength2, 2));
pres3 = 612*ones(1,size(InflatedLength3, 2));
% pres = [pres1 pres3 pres4];
pres = [pres1 pres2 pres3 pres4];

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
    
%     TorqueHand1(:,:,i) = TorqueHand(1,1:size(TorqueHand1,2),i);
%     TorqueHand3(:,:,i) = TorqueHand(1,(size(TorqueHand1,2)+1):(size(TorqueHand1,2)+size(TorqueHand3,2)),i);
%     TorqueHand4(:,:,i) = TorqueHand(1,((size(TorqueHand1,2)+size(TorqueHand3,2))+1):size(TorqueHand,2),i);
    TorqueHand1(:,:,i) = TorqueHand(1,1:size(TorqueHand1,2),i);
    TorqueHand2(:,:,i) = TorqueHand(1,(size(TorqueHand1,2)+1):(size(TorqueHand1,2)+size(TorqueHand2,2)),i);
    TorqueHand3(:,:,i) = TorqueHand(1,(size(TorqueHand1,2)+size(TorqueHand2,2)+1):(size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2)),i);
    TorqueHand4(:,:,i) = TorqueHand(1,((size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2))+1):size(TorqueHand,2),i);

end

%% Prepare for plotting

% % Since one of the hand measurement recordings was corrupted, use a different 'Angle' for tests 1+3+4
% AngleX = [Angle1 Angle3 Angle4];
% TorqueHandX = [TorqueHand1 TorqueHand3 TorqueHand4];

% Since two of the hand measurement recordings was corrupted, use a different 'Angle' for tests 3+4
AngleX = [Angle3 Angle4];
TorqueX = [Torque3 Torque4];
TorqueHandX = [TorqueHand3 TorqueHand4];

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
Ma = Vas_Pam_48cm.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
ax11 = subplot(2,1,1);
hold on
title('\bf Expected vs measured $r_{\hat{k}}$', 'Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf Z axis $r_{\hat{k}}$, m','Interpreter','latex')
pp21 = plot(phiD,G,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected $r_{\hat{k}}$');
ss2 = scatter(AngleX, ICRtoMuscleX,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured $r_{\hat{k}}$');
% ss2 = scatter(Angle, ICRtoMuscle,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured $r_{\hat{k}}$');
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
ss2_1 = scatter(Angle1, ICRtoMuscle1,sz,'filled','MarkerFaceColor',c1,'DisplayName','\bf Measured $r_{\hat{k}}$, JM LC');
ss2_2 = scatter(Angle2, ICRtoMuscle2,sz,'filled','MarkerFaceColor',c2,'DisplayName','\bf Measured $r_{\hat{k}}$, JM LC');
ss2_3 = scatter(Angle3, ICRtoMuscle3,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured $r_{\hat{k}}$, BB FS');
ss2_4 = scatter(Angle4, ICRtoMuscle4,sz,'filled','MarkerFaceColor',c6,'DisplayName','\bf Measured $r_{\hat{k}}$, BB LC');
set(ax12,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial')
lgdMa12 = legend('Interpreter','latex');
lgdMa12.FontSize = 12;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values
strain = (rest-(Vas_Pam_48cm.MuscleLength-tendon-2*fitting))/rest;
relstrain = (strain)./KMAX;
realRel = (rest-InflatedLength)/rest/KMAX;
realRelX = (rest-InflatedLengthX)/rest/KMAX;
realRel1 = (rest-InflatedLength1)/rest/KMAX;
realRel2 = (rest-InflatedLength2)/rest/KMAX;
realRel3 = (rest-InflatedLength3)/rest/KMAX;
realRel4 = (rest-InflatedLength4)/rest/KMAX;

figure
ax21 = subplot(2,1,1);
hold on
title('Expected vs measured \epsilon^*','Interpreter','tex')
xlabel('Knee angle, \circ','Interpreter','tex')
ylabel('\epsilon^*','Interpreter','tex')
plot(phiD,relstrain,'Color',c7,'LineWidth',2,'DisplayName','\bf Expected \epsilon^*')
% scatter(AngleX,realRel,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured \epsilon^*')
scatter(AngleX,realRelX,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured \epsilon^*')
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
scatter(Angle1,realRel1,sz,'filled','MarkerFaceColor',c1,'DisplayName','\bf Measured \epsilon^*, JM LC')
scatter(Angle2,realRel2,sz,'filled','MarkerFaceColor',c2,'DisplayName','\bf Measured \epsilon^*, JM LC')
scatter(Angle3,realRel3,sz,'filled','MarkerFaceColor',c4,'DisplayName','\bf Measured \epsilon^*, BB FS')
scatter(Angle4,realRel4,sz,'filled','MarkerFaceColor',c6,'DisplayName','\bf Measured \epsilon^*, BB LC')
set(ax22,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdEp1 = legend('Interpreter','tex');
lgdEp1.FontSize = 12;
hold off

%% Plot measured versus expected BPA length
MuscleLength = Vas_Pam_48cm.MuscleLength-2*fitting-tendon;

figure
hold on
title('\bf Expected vs measured $l_{m}$','Interpreter','latex')
xlabel('\bf Knee angle, \circ','Interpreter','tex')
ylabel('\bf $l_{m}$, m','Interpreter','latex')
plot(phiD,MuscleLength,'DisplayName','\bf Expected $l_{m}$')
scatter(AngleX,InflatedLengthX,'DisplayName','\bf Measured $l_{m}$')
% scatter(Angle,InflatedLength,'DisplayName','\bf Measured $l_{m}$')
set(ax22,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2, 'FontName','Arial')
lgdLm = legend('Interpreter','latex');
lgdLm.FontSize = 12;
hold off

%% Plotting the Torque results
figure
hold on
title('Iso. Torque vs {\theta_{k}}, Pinned, Extensor, l_{rest} = 48 cm','interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','interpreter','tex')
gca1 = gca;
gcf1 = gcf;
set(gca1,'FontSize', 12, 'FontWeight', 'bold','LineWidth',2,'FontName','Arial','XLim',[-125 40])

PL1 = plot(phiD, Theoretical,'Color',c4,'Linewidth',2,'DisplayName','Theoretical');
scM = scatter(Angle,Torque,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Measured');
scH = scatter(AngleX,TorqueHandX(:,:,4),sz,'filled','MarkerFaceColor',c1,'DisplayName','Back calculated');
% scH = scatter(Angle,TorqueHand(:,:,4),sz2,'filled','MarkerFaceColor',c1,'DisplayName','Back calculated');

lgd1 = legend;
hold off

%% Plotting Torque tests, Color results by tests

figure
hold on
title('Iso. Torque vs {\theta_{k}}, Pinned, Extensor, l_{rest} = 48 cm','interpreter','tex')
xlabel('Knee angle, \circ','FontWeight','bold','interpreter','tex')
ylabel('Torque, N{\cdot}m','FontWeight','bold','interpreter','tex')
gca2 = gca;
gcf2 = gcf;
set(gca2,'FontSize', 12, 'FontWeight','bold','LineWidth',2,'FontName','Arial','XLim',[-125 40])

PL2 = plot(phiD, Theoretical,'Color',c4,'Linewidth',2,'DisplayName','Theoretical');

scM1 = scatter(Angle1,Torque1,sz,'d','filled','MarkerFaceColor',c1,'DisplayName','JM LC');
scM2 = scatter(Angle2,Torque2,sz,'d','filled','MarkerFaceColor',c3,'DisplayName','JM FS');
scM3 = scatter(Angle3,Torque3,sz,'d','filled','MarkerFaceColor',c5,'DisplayName','BB FS');
scM4 = scatter(Angle4,Torque4,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','BB LC');
scH1 = scatter(Angle1,TorqueHand1(:,:,4),sz2,'filled','MarkerFaceColor',c1,'DisplayName','JM LC hand');
scH2 = scatter(Angle2,TorqueHand2(:,:,4),sz2,'filled','MarkerFaceColor',c3,'DisplayName','JM FS hand');
scH3 = scatter(Angle3,TorqueHand3(:,:,4),sz2,'filled','MarkerFaceColor',c5,'DisplayName','BB FS hand');
scH4 = scatter(Angle4,TorqueHand4(:,:,4),sz2,'filled','MarkerFaceColor',c7,'DisplayName','BB LC hand');

lgd2 = legend;
hold off

%% Mean and RMSE
Tqz = cell(1,1);
Tqz{1} = Vas_Pam_48cm.Torque(:,3,4);    	%Calculated Torque, new simplified exponential equation w/o optimized fitting length
%Tqz{2} = Vas_Pam_48cm_adj.Torque(:,3,4);   %Calculated Torque, adjusted with optimized fitting length
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
