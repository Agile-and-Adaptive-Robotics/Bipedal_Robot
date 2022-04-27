%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; %close all;

f = fullfile('github/bipedal_robot/code/matlab');
qt = addpath(genpath(f));

restingLength = 0.415; %resting length, m
kmax = 0.350; %Length at maximum contraction, m
dia = 10;

load KneeFlx_10mm_42cm.mat
Theoretical = TorqueR(:,3)';
%rest = 0.415, tendon = 0.012

%% Adjustment to Theoretical Calculation
rest = 0.5261;
tendon = 0.0125;
kmax = 0.3506;
Bifemsh_Pam_adj = MonoPamDataExplicit(Name, Location, CrossPoint, dia, T_Pam, rest, kmax, tendon, fitting, pres);
Theo_adj = Bifemsh_Pam_adj.Torque(:,3)';

rest = 0.415; %set resting length back to measured value
tendon = 0.012; %set tendon back to measured(?) value
kmax = 0.350; %set kmax back to measured value
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

F = zeros(1,size(InflatedLength, 2));

TorqueHand = zeros(1,size(InflatedLength, 2));

%load pressure where applicable
runsperseries = 10;

pres = zeros(1,runsperseries);

     for j = 1:runsperseries(i)
                file_name = sprintf('FlxTest%0.0f.mat', j);
                load(file_name,'Stats')
                pres(1,j) = Stats{'Mean',2};
     end

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, dia, pres(i), kmax);
    TorqueHand(i) = -ICRtoMuscle(i)*F(i);  %Moment will be negative because it causes flexion
end

%% Mean and RMSE
X1 = linspace(min(Angle),max(Angle),size(Angle,2));      %Range of motion
mod = 'poly4';
fitOptions = fitoptions(mod, 'Normalize', 'on','Robust','on');
[mdl1u, gof1] = fit(Angle',Torque',mod,fitOptions)
TorqueStdu = gof1.rmse;
TorqueMeanu = feval(mdl1u,X1)';


modp = 'poly3';
fitOp = fitoptions(modp,'Normalize','on','Robust','on');
[mdl1, gofp1] = fit(Angle',Torque',modp,fitOp)
TorqueStd = gofp1.rmse;
TorqueMean = feval(mdl1,X1)';

[mdl2, gofp2] = fit(Angle',TorqueHand',modp,fitOp);
HandStd = gofp2.rmse;
HandMean = feval(mdl2,X1)';

%% Plotting with polynomial solver
%Matlab hex color values:
c1 = '#FFD700'; %gold
c2 = '#FFB14E'; %orange
c3 = '#FA8775'; %light orange
c4 = '#EA5F94'; %pink
c5 = '#CD34B5'; %magenta
c6 = '#9D02D7'; %magenta 2
c7 = '#0000FF'; %indigo
sz = 60;        %size of data points

%close(get(gcf,'Number'));
figure('units','normalized','position',[0.0892 0.1100 0.5317 0.8150])
hold on
title('Isometric Torque for Biomimetic Knee, 10mm Flexor, 41.5cm long')
xlabel('Knee angle, degrees, Flexion(-),Extension(+)')
ylabel('Torque, N*m')
gca1 = gca;
gcf1 = gcf;
% set(gcf,'Position',[1 384 950 612]);
% set(gca,'FontSize', 12, 'FontWeight', 'bold','XMinorGrid','on','XMinorTick','on','YMinorGrid','on','YMinorTick','on');
set(gca,'FontSize', 12, 'FontWeight', 'bold');
plot(phiD, Theoretical,'Color',c5,'Linewidth',2,'DisplayName','Expected')
chr = ['|Adjustment, mm' newline '|rest-5' newline '|tendon+26' newline '|kmax-4'];
plot(phiD, Theo_adj,'Color',c3,'Linewidth',2,'DisplayName',chr)

% Xnew=[X,fliplr(X)];
% Y1=[TorqueMean+TorqueStd,fliplr(TorqueMean-TorqueStd)];
% Y2=[HandMean+HandStd,fliplr(HandMean-HandStd)];
% fill(Xnew,Y1,[1 0.4 0.8],'DisplayName','Scale torque SD','FaceAlpha',0.25);
% fill(Xnew,Y2,[.6 1.0 .6],'DisplayName','Hand torque SD','FaceAlpha',0.25);

plot(X1,TorqueMean,'--','Color',c7,'Linewidth',2,'DisplayName','Torque mean, scale')
plot(X1,HandMean,'--','Color',c1,'Linewidth',2,'DisplayName','Torque expected, hand measurements')

sc1 = scatter(Angle,Torque,sz,'d','filled','MarkerFaceColor',c7,'DisplayName','Torque data, experimental');
sc2 = scatter(Angle,TorqueHand,sz,'filled','MarkerFaceColor',c1,'DisplayName','Torque data, Hand calc');

% Legend = legend('Location','eastOutside');
legend
hold off

%% Plot error and standard deviation as bar graphs
% xb=categorical({'Theoretical','Scale','Hand'});
% xb = reordercats(xb,{'Theoretical','Scale','Hand'});
% yb = [mean(Theoretical) mean(TorqueMean) mean(TorqueHand)];
% std_dev = [0 mean(TorqueStd) mean(HandStd)];
% figure
% hold on
% b = bar(xb,yb,'FaceColor',[0 0.8 1.0]);
% gca3 = gca;
% gcf3 = gcf;
% ylabel('Torque, N*m')
% title('Mean Torque Values and SD for Each Calculation Method')
% set(gcf,'Position',[0 0 950 612]);
% set(gca,'FontSize', 18, 'FontWeight', 'bold');
% b.CData = [0 0.4470 0.7410; 0  0  0; 1  0  0];
% errb = errorbar(yb,std_dev ,'LineStyle','none','LineWidth',4,'CapSize',20);
% hold off