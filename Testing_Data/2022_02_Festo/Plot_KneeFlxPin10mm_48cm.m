%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.485; %resting length, m
kmax = 0.398; %Length at maximum contraction, m

load KneeFlxPin_10mm_48cm.mat
Theoretical = cell(2, size(TorqueR,3));
Theoretical{1,1} = 'Hunt Eq.';
Theoretical{1,2} = 'Exponential Eq.';
Theoretical{1,3} = 'Polynomial Eq.';
Theoretical{1,4} = 'Exponential Eq., Simplified';
Theoretical{1,5} = 'Polynomial Eq., Simplified';
for i = 1:length(Theoretical)
    Theoretical{2,i} = TorqueR(:,3,i)';
end
%% Test 1&2 done with CALT load cell. Tests 3-5 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 
%Test 1 == sheet FlxTest10mm_1 from Results_table10mm_pinned_LoadCell
%Test 2 == sheet FlxTest10mm_2 from Results_table10mm_pinned_LoadCell
%Test 3 == sheet FlxTest10mm_1 from Results_table10mm_pinned_FishScale
%Test 4 == sheet FlxTest10mm_4 from Results_table10mm_pinned_FishScale
%Test 5 == sheet FlxTest10mm_2 from Results_table10mm_pinned_FishScale [not
%%used due to uncertainty about resting length]
%% Torque calculated from measurements

Angle1 = [-5.5	-10.5	-25.5	-4.5	-14.5	-20	-33	-39.5	-43	-50	-62.5	-60	-55	-53.5	-46	-40	-36.5	-33.5	-29.5	-24.5	-17.5	-14.5	-5.5	-3	-0.5	3.5	11.5	19.5];
Angle2 = [6	5.5	-1	-15	-26.6	-35	-42.5	-56	-56	-63	-64	-55	-48	-42	-35	-26	-19	-11	0.5	11.5];
Angle3 = [-10	-14	-24	-27	-23	-21	-35	-45	-52	-61	-65];
Angle4 = [-22	-25	-32	-35	-37	-42	-44.5	-52	-57	-62	-68	-71	-69.5	-67	-61.5	-48	-33	-34	-25.5	-36	-2	-8	-13.5	-19.5	-24	-24	-20	-19	-15	-13	-2];
%Angle5 = [-18.5	-30.5	-39.5	-36.5	-42	-46.5	-52.5	-69.5	-68.5	-69.9	-60.5	-53	-47.5	-37.5	-44.5	-31.5	-26.5	-23];
Angle = [Angle1, Angle2, Angle3, Angle4];

Torque1 = [-13.55019944	-12.71574633	-12.2274503	-13.57429518	-12.61783318	-11.4271733	-7.828935983	-6.496239587	-4.940487337	-2.826134001	-1.210800042	-0.110215947	-1.513744872	-3.192740824	-4.861595823	-6.86945259	-8.599788984	-10.3582314	-11.8636379	-13.47755504	-13.77813076	-14.88459711	-14.53453947	-14.33814406	-14.02122399	-12.53904665	-11.16785477	-8.28360707];
Torque2 = [-12.02448648	-12.01072043	-12.46707433	-11.42926509	-9.361208447	-7.676179585	-4.957997672	-1.67225365	-0.807464392	-0.003574341	-2.019875	-4.170597119	-6.466679019	-8.166509485	-10.60330251	-12.64002428	-13.93045057	-14.68948292	-14.67482254	-13.15130291];
Torque3 = [-13.93505973	-13.41723648	-11.54052398	-9.443681562	-13.61674065	-13.66679506	-8.935879493	-5.788289674	-3.353532947	-1.833248329	-0.531376327];
Torque4 = [-15.24648357	-13.94170519	-12.36216232	-11.04703867	-10.0579279	-8.249812691	-6.962421539	-4.874859191	-3.096568251	-2.046884098	-1.0207759	0	0	-1.0207759	-2.834147213	-7.995951798	-11.5993906	-13.3993773	-15.9919036	-11.86241533	-16.23568514	-17.29899289	-17.10478208	-16.01872483	-14.22288584	-14.88045526	-15.17480498	-16.23203091	-16.25073881	-17.09230224	-16.15641241];
%Torque5 = [-19.06895474	-15.75701078	-16.75788099	-14.15214602	-10.75456087	-8.935879493	-5.252336926	-3.397585151	-1.308275082	-1.5822075	-3.390835418	-6.554396551	-8.661915638	-13.68558121	-16.00098846	-17.05641248	-18.14010372	-18.66468509];
Torque = [Torque1, Torque2, Torque3, Torque4];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength1 = [461.5	456	450.5	460.5	452	449	432	428	423	415	409.5	406	410	414.5	420.5	426	430	433.5	441	445	451	454	455.5	458.5	462	461.5	465	467]/1000;
InflatedLength2 = [467	466	463	450	443	431	423	411	408	404	407	413	420	422	433	442	448	455	460	466]/1000;
InflatedLength3 = [460	457	451	445	450	452	438	427	420	415	409]/1000;
InflatedLength4 = [450	445	439	435	431	426	421	410	408	406	402	398	398	400	404	420	429	434	444	431	461	459	453	450	442	442	443	445	449	454	458]/1000;
%InflatedLength5 = [459	448	449	433	423	419.5	410.5	406	396	392.5	401.5	407	411.5	426	432	436	441	442.5]/1000;
InflatedLength = [InflatedLength1, InflatedLength2, InflatedLength3, InflatedLength4];

ICRtoMuscle1 = [62	74	76	64	72	76	79.5	81.5	77	81	74	68	74	77	75	79	80	78	79	77	75.5	72	61.5	62.5	54.5	48	38.5	30]/1000;
ICRtoMuscle2 = [45	50	55	65	70	74	73	72	65	63	66	71	75	75	76	75	69	65	56	44]/1000;
ICRtoMuscle3 = [66	70	75	77	74	75	78	77	75	68	69]/1000;
ICRtoMuscle4 = [75	77	80	83	84	84	84	79	77	75	69	65	65	67	73	84	83	82	79	83	58	64	69	73	78	78	75	75	75	70	62]/1000;
%ICRtoMuscle5 = [76	84	81	85	84	83	81	85	64	66.5	70	81	83.5	85.5	85	86	80.5	79.5]/1000;
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle2, ICRtoMuscle3, ICRtoMuscle4];

F1 = zeros(1,size(InflatedLength1, 2));
F2 = zeros(1,size(InflatedLength2, 2));
F3 = zeros(1,size(InflatedLength3, 2));
F4 = zeros(1,size(InflatedLength4, 2));
%F5 = zeros(1,size(InflatedLength5, 2));
F = [F1, F2, F3, F4];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand3 = zeros(1,size(InflatedLength3, 2));
TorqueHand4 = zeros(1,size(InflatedLength4, 2));
%TorqueHand5 = zeros(1,size(InflatedLength5, 2));
TorqueHand = [TorqueHand1, TorqueHand2, TorqueHand3, TorqueHand4];

%load pressure where applicable
test = [1 2];
runsperseries = [28 19];

pres1 = zeros(1,runsperseries(1));
pres2 = zeros(1,runsperseries(2));
 
 for i = 1:size(test, 2)   
   if i == 1 
     for j = 1:runsperseries(i)
                file_name = sprintf('FlexTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2};
     end
   else  
    for j = 1:runsperseries(i)
                file_name = sprintf('FlexTest%0.0f_%0.0f.mat', test(i),j);
                load(file_name,'Stats')
                pres2(1,j) = Stats{'Mean',2};
    end
   end
 end
presx = mean(pres2);                             %Missing 20th data point for Test 2 
pres3 = 612*ones(1,size(InflatedLength3, 2));
pres4 = 612*ones(1,size(InflatedLength4, 2));
%pres5 = 612*ones(1,size(InflatedLength5, 2));
pres = [pres1 pres2 presx pres3 pres4];

for i = 1:size(InflatedLength, 2)  
    F(1,i,1) = festo3(InflatedLength(i), restingLength, 10, pres(i), kmax);    
    TorqueHand(i) = -ICRtoMuscle(i)*F(i);  %Torque will be negative because it is causing flexion
end

KMAX = (restingLength-kmax)/restingLength;
rel = ((restingLength-InflatedLength)/restingLength)/KMAX;
Fn = bpaForce10(restingLength,rel,pres);

for i = 2:(size(Fn,3)+1)
    F(:,:,i) = Fn(:,:,i-1);
end

for i = 1:size(F,3)
    TorqueHand(:,:,i) = -ICRtoMuscle.*F(1,:,i);  %Torque will be negative because it is causing flexion
    
    TorqueHand1(:,:,i) = TorqueHand(1,1:size(TorqueHand1,2),i);
    TorqueHand2(:,:,i) = TorqueHand(1,(length(TorqueHand1)+1):(size(TorqueHand1,2)+size(TorqueHand2,2)),i);
    TorqueHand3(:,:,i) = TorqueHand(1,(size(TorqueHand1,2)+size(TorqueHand2,2)+1):(size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2)),i);
    TorqueHand4(:,:,i) = TorqueHand(1,(size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2)+1):(size(TorqueHand1,2)+size(TorqueHand2,2)+size(TorqueHand3,2)+size(TorqueHand4,2)),i);
end


%% Process data for further curve fitting

Hand1 = TorqueHand(:,:,1); %'Hunt Eq.';
Hand2 = TorqueHand(:,:,2); %'Exponential Eq.';
Hand3 = TorqueHand(:,:,3); %'Polynomial Eq.';
Hand4 = TorqueHand(:,:,4); %'Exponential Eq., Simplified';
Hand5 = TorqueHand(:,:,5); %'Polynomial Eq., Simplified';

A = [Angle', Torque'];                              % Arrange Data
[UA,~,idx] = unique(A(:,1));                        % Find repeated angle values
NEW_A = [UA,accumarray(idx,A(:,2),[],@mean)];       % Take mean of torque at repeated angle values

NEW_ang = NEW_A(:,1);           % New angle
NEW_tq = NEW_A(:,2);            % New torque

SSE_1 = sum((Torque-Hand1).^2,'omitnan')
SSE_2 = sum((Torque-Hand2).^2,'omitnan')
SSE_3 = sum((Torque-Hand3).^2,'omitnan')
SSE_4 = sum((Torque-Hand4).^2,'omitnan')
SSE_5 = sum((Torque-Hand5).^2,'omitnan')
RMSE_1 = sqrt(sum((Torque-Hand1).^2,'omitnan')/length(Torque))
RMSE_2 = sqrt(sum((Torque-Hand2).^2,'omitnan')/length(Torque))
RMSE_3 = sqrt(sum((Torque-Hand3).^2,'omitnan')/length(Torque))
RMSE_4 = sqrt(sum((Torque-Hand4).^2,'omitnan')/length(Torque))
RMSE_5 = sqrt(sum((Torque-Hand5).^2,'omitnan')/length(Torque))

%% Plot expected versus measured moment arm
Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque

figure
hold on
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
title('Expected vs measured moment arm')
xlabel('Knee angle, degrees')
ylabel('Moment Arm, z axis (m)')
pp = plot(phiD,G,'DisplayName','MA expected');
ss = scatter(Angle, ICRtoMuscle,'DisplayName','MA measured');
lgdMa = legend;
lgdMa.FontSize = 10;
hold off

figure
hold on
fig1 = gcf;
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
title('Expected vs measured moment arm')
xlabel('Knee angle, \circ')
ylabel('Moment Arm, z axis (m)')

subplot(2,2,1)
hold on
pp1 = plot(phiD,G,'DisplayName','MA expected');
ss1 = scatter(Angle1, ICRtoMuscle1,'DisplayName','MA measured, JM LC');
title('Expected vs measured moment arm')
xlabel('Knee angle, \circ')
ylabel('Moment Arm, z axis (m)')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 10;
hold off

subplot(2,2,2)
hold on
pp2 = plot(phiD,G,'DisplayName','MA expected');
ss2 = scatter(Angle2, ICRtoMuscle2,'DisplayName','MA measured, BB LC');
title('Expected vs measured moment arm')
xlabel('Knee angle, \circ')
ylabel('Moment Arm, z axis (m)')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 10;
hold off

subplot(2,2,3)
hold on
pp3 = plot(phiD,G,'DisplayName','MA expected');
ss3 = scatter(Angle3, ICRtoMuscle3,'DisplayName','MA measured, BB FS');
title('Expected vs measured moment arm')
xlabel('Knee angle, \circ')
ylabel('Moment Arm, z axis (m)')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 10;
hold off

subplot(2,2,4)
hold on
pp4 = plot(phiD,G,'DisplayName','MA expected');
ss4 = scatter(Angle4, ICRtoMuscle4,'DisplayName','MA measured, BB FS');
title('Expected vs measured moment arm')
xlabel('Knee angle, \circ')
ylabel('Moment Arm, z axis (m)')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 10;
hold off

%% Plot relative strain versus angle. Compare strain, relative strain, and measured values

% tendon = 0;
% fitting = 0.40;
strain = (restingLength-(Bifemsh_Pam.MuscleLength-tendon-2*fitting))/restingLength;
relstrain = (strain)./KMAX;
realRel = (restingLength-InflatedLength)/restingLength/KMAX;

figure
hold on
plot(phiD,relstrain,'DisplayName','Expected Relative Strain')
scatter(Angle,realRel,'DisplayName','Measured Relative Strain')
title('Expected vs measured relative strain')
xlabel('Knee angle, \circ')
ylabel('strain/kmax')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 10;
hold off

%% Plot measured versus expected strain (like above, but not normalized
realStrain = (restingLength-InflatedLength)/restingLength;

figure
hold on
plot(phiD,strain,'DisplayName','Expected Strain')
scatter(Angle,realStrain,'DisplayName','MeasuredStrain')
title('Expected vs measured strain')
xlabel('Knee angle, \circ')
ylabel('strain')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 10;
hold off

%% Plot measured versus expected strain (like above, but not normalized
MuscleLength = Bifemsh_Pam.MuscleLength-2*fitting-tendon;

figure
hold on
plot(phiD,MuscleLength,'DisplayName','Expected Muscle Length')
scatter(Angle,InflatedLength,'DisplayName','Measured Length')
title('Expected vs measured muscle length')
xlabel('Knee angle, \circ')
ylabel('Length, m')
ax1 = gca;
ax1.FontSize = 12;
ax1.FontWeight = 'bold';
ax1.FontName = 'Arial';
ax1.YAxis.LineWidth = 2; ax1.YAxis.FontSize = 10;
ax1.XAxis.LineWidth = 2; ax1.XAxis.FontSize = 10;
lgdMa = legend;
lgdMa.FontSize = 10;
hold off

%% Plotting measured versus expected Torque values solver
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

PL = cell(1, size(Theoretical,2));
sc = cell(1, size(Theoretical,2));
sc1 = cell(1, size(Theoretical,2));
sc2 = cell(1, size(Theoretical,2));
sc3 = cell(1, size(Theoretical,2));
sc4 = cell(1, size(Theoretical,2));
Disp = cell(1, 2*size(Theoretical,2));

figure

scM = scatter(Angle,Torque,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Torque data, measured');
hold on
    for i = 1:size(Theoretical,2)
        txt = Theoretical{1,i};
        T1 = 2*i-1;
        H1 = 2*i;
        Disp{T1} = sprintf('Theoretical Torque, %s',txt);
        Disp{H1} = sprintf('A posteriori torque, %s',txt);
        PL{i} = plot(phiD, Theoretical{2,i},'Color',c{i},'Linewidth',2,'DisplayName',Disp{T1});
        sc{i} = scatter(Angle,TorqueHand(:,:,i),sz,'filled','MarkerFaceColor',c{i},'DisplayName',Disp{H1});
    end
title('Isometric Torque vs Knee Angle, 10mm Flexor, 48.5cm long','FontSize',12,'FontWeight','Bold')
xlabel('Knee Flexion(-)/Extension(+), \circ','FontSize',12)
ylabel('Torque, N \bullet m','FontSize',12)
ax2 = gca;
ax2.FontSize = 12;
ax2.FontWeight = 'bold';
ax2.FontName = 'Arial';
ax2.YAxis.LineWidth = 2; ax2.YAxis.FontSize = 10;
ax2.XAxis.LineWidth = 2; ax2.XAxis.FontSize = 10;
lgd = legend;
lgd.FontSize = 10;
lgd.Location = 'southwest';
hold off

figure
t = tiledlayout(2,2);
t.Title.String = 'Isometric Torque vs Knee Angle, 10mm Flexor, 48.5cm long';
t.Title.FontSize =12;
t.Title.FontWeight = 'bold';
t.Title.FontName = 'Arial';
t.XLabel.String = 'Knee Flexion(-)/Extension(+), \circ';
t.XLabel.FontSize =12;
t.XLabel.FontWeight = 'bold';
t.XLabel.FontName = 'Arial';
t.YLabel.String = 'Torque, N \bullet m';
t.YLabel.FontSize =12;
t.YLabel.FontWeight = 'bold';
t.YLabel.FontName = 'Arial';


ax4 = nexttile;
scM1 = scatter(Angle1,Torque1,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Torque data, measured');
hold on
    for i = 1:size(Theoretical,2)
        txt = Theoretical{1,i};
        T1 = 2*i-1;
        H1 = 2*i;
        Disp{T1} = sprintf('Theoretical Torque, %s',txt);
        Disp{H1} = sprintf('A posteriori torque, %s',txt);
        PL{i} = plot(phiD, Theoretical{2,i},'Color',c{i},'Linewidth',2,'DisplayName',Disp{T1});
        sc1{i} = scatter(Angle1,TorqueHand1(:,:,i),sz,'filled','MarkerFaceColor',c{i},'DisplayName',Disp{H1});
    end
title('Subplot 1: JM LC/hand')
ax4.FontSize = 10;
ax4.FontWeight = 'bold';
ax4.FontName = 'Arial';
ax4.YAxis.LineWidth = 2; ax4.YAxis.FontSize = 10;
ax4.XAxis.LineWidth = 2; ax4.XAxis.FontSize = 10;
lgd1 = legend;
lgd1.FontWeight = 'bold';
lgd1.FontSize = 8;
lgd1.FontName = 'Arial';
hold off

ax5 = nexttile;
scM2 = scatter(Angle2,Torque2,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Torque data, measured');
hold on
    for i = 1:size(Theoretical,2)
        txt = Theoretical{1,i};
        T1 = 2*i-1;
        H1 = 2*i;
        Disp{T1} = sprintf('Theoretical Torque, %s',txt);
        Disp{H1} = sprintf('A posteriori torque, %s',txt);
        PL{i} = plot(phiD, Theoretical{2,i},'Color',c{i},'Linewidth',2,'DisplayName',Disp{T1});
        sc2{i} = scatter(Angle2,TorqueHand2(:,:,i),sz,'filled','MarkerFaceColor',c{i},'DisplayName',Disp{H1});
    end

title('Subplot 2: BB LC/hand')
ax5.FontSize = 10;
ax5.FontWeight = 'bold';
ax5.FontName = 'Arial';
ax5.YAxis.LineWidth = 2; ax5.YAxis.FontSize = 10;
ax5.XAxis.LineWidth = 2; ax5.XAxis.FontSize = 10;
lgd2 = legend;
lgd2.FontSize = 8;
lgd2.FontWeight = 'bold';
hold off

ax6 = nexttile;
scM3 = scatter(Angle3,Torque3,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Torque data, measured');
hold on
    for i = 1:size(Theoretical,2)
        txt = Theoretical{1,i};
        T1 = 2*i-1;
        H1 = 2*i;
        Disp{T1} = sprintf('Theoretical Torque, %s',txt);
        Disp{H1} = sprintf('A posteriori torque, %s',txt);
        PL{i} = plot(phiD, Theoretical{2,i},'Color',c{i},'Linewidth',2,'DisplayName',Disp{T1});
        sc3{i} = scatter(Angle3,TorqueHand3(:,:,i),sz,'filled','MarkerFaceColor',c{i},'DisplayName',Disp{H1});
    end

title('Subplot 3: BB FS/hand')
ax6.FontSize = 12;
ax6.FontWeight = 'bold';
ax6.FontName = 'Arial';
ax6.YAxis.LineWidth = 2; ax6.YAxis.FontSize = 10;
ax6.XAxis.LineWidth = 2; ax6.XAxis.FontSize = 10;
lgd3 = legend;
lgd3.FontSize = 8;
lgd3.FontWeight = 'bold';
hold off

ax7 = nexttile;
scM4 = scatter(Angle4,Torque4,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Torque data, measured');
hold on
    for i = 1:size(Theoretical,2)
        txt = Theoretical{1,i};
        T1 = 2*i-1;
        H1 = 2*i;
        Disp{T1} = sprintf('Theoretical Torque, %s',txt);
        Disp{H1} = sprintf('A posteriori torque, %s',txt);
        PL{i} = plot(phiD, Theoretical{2,i},'Color',c{i},'Linewidth',2,'DisplayName',Disp{T1});
        sc4{i} = scatter(Angle4,TorqueHand4(:,:,i),sz,'filled','MarkerFaceColor',c{i},'DisplayName',Disp{H1});
    end
title('Subplot 4: BB FS/hand')
ax7.FontSize = 12;
ax7.FontWeight = 'bold';
ax7.FontName = 'Arial';
ax7.YAxis.LineWidth = 2; ax7.YAxis.FontSize = 10;
ax7.XAxis.LineWidth = 2; ax7.XAxis.FontSize = 10;
lgd4 = legend;
lgd4.FontSize = 8;
lgd4.FontWeight = 'bold';
hold off


