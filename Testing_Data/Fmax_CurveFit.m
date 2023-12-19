%% Look for equations of max force and max strain as functions of resting length
%  Data from experiments. See Excel spreadsheet "StraightForceTest"

clear; clc; close all

A = [845	840     785     780     710     709        571     112     415     453     490     518	551     361     54      27      69      275     151     193     10      120     220     260     281     281     83      132     182     233     281     382     426     521;
    0.1716	0.1726	0.1694	0.1705	0.1775	0.1749	0.1821	0.1607	0.1663	0.1611	0.1878	0.1660	0.1742	0.1717	0.1481	0.1481	0.1594	0.1709	0.1523	0.1658	0.1000	0.1667	0.1591	0.1615	0.1637	0.1459	0.1566	0.1591	0.1648	0.1717	0.1744	0.1780	0.1761	0.1612;
    447.1	447.1	472.0	472.0	452.3	452.3	461.6	334.74	445        460	458.6	455.83	453.1	436.4	238.2	135.3	271.5	412.3	347.2	396.2	8.3     341.27	377.75	383.53	419.44	407.13	317.15	363.61	380.25	394.95	436.00	459.11	455.37	484.23];


restingL = A(1,1:23)/1000;  %Resting length

kmax = A(2,1:23);          %Maximum strain

Fmax = A(3,1:23);          %Maximum force
%Fmax_norm = Fmax/max(Fmax);

%% Nonlinear fit for Fmax
modelfun = @(b,xm)b(1)*620*atan(b(2)*620*(xm-0.0075)); %based on the shape, it looks like resting length reaches a limit
%modelfun = @(b,xm)-exp(-b(1)*(xm-b(2))+b(3))+b(4);
beta0 = [0.4909 0.0317];
%beta0 = [11 0.41 1.2 478];
opts = statset('fitnlm');
opts.MaxIter = 1000;
opts.Display = 'final';
opts.RobustWgtFun = 'logistic';
mdl = fitnlm(restingL,Fmax,modelfun,beta0,'Exclude',[],'Options',opts)
x = linspace(0,1);
p1 = feval(mdl,x);

md2 = @(x) 0.4909*620*atan(0.0317*620*(x-0.0075));     %from Fmax_fitTool
p2 = feval(md2,x);

figure
plot(restingL,Fmax,'o',x,p2,'--')
title('F_{max10} vs. l_{rest}')
xlabel('l_{rest}, m')
ylabel('F_{max}, N')
legend('Data','Model Fit')

%% Fit maximum force (i.e. zero contraction BPA) for different resting lengths and pressures
%Load data from StraightForceTest Fmax(L) tab. It uses data from some of
%Lawrence's tests.

% Force (z-axis) unsorted
Zun = [840	0	NaN	114.32	191.8	275.2	355.4	447.1;
    780	0	NaN	139.2	216	299.7	378.1	472;
    709	0	NaN	124.21	201.4	276.02	355.9	452.32;
    571	0	NaN	136.1	212.39	288.2	366.94	461.6;
    112	0	NaN	NaN	NaN	NaN	NaN	334.74;
    415	0	NaN	NaN	NaN	NaN	NaN	445;
    455	0	NaN	NaN	NaN	NaN	NaN	460;
    490	0	NaN	NaN	NaN	NaN	NaN	458.6;
    518	0	NaN	NaN	NaN	NaN	NaN	455.83;
    551	0	NaN	129.26	201.3	282.34	360.24	453.14;
    361	0	NaN	121.5	192.91	269	344.5	436.4;
    54	0	NaN	56.1	103.37	155.45	207	238.2;
    27	0	NaN	NaN	NaN	NaN	NaN	135.32;
    69	0	NaN	NaN	NaN	NaN	NaN	271.48;
    275	0	NaN	NaN	NaN	NaN	NaN	412.3;
    151	0	NaN	NaN	NaN	NaN	NaN	347.11;
    193	0	NaN	NaN	NaN	NaN	NaN	396.17;
    10	0	NaN	NaN	NaN	NaN	NaN	8.3;
    120	0	11.2	71.6	132.0	194.3	260.6	341.27;
    220	0	19.9	83.8	147.8	217.4	289.2	377.75;
    260	0	30.6	92.7	156.1	225.1	296.1	383.53;
    281	0	54.4	120.7	187.5	258.1	330.3	419.44;
    281	0	37.7	104.7	173.0	244.5	317.6	407.13;
%     83	0	18.2	67.2	120.7	180.3	242.2	317.2;
%     132	0	28.6	83.0	144.1	210.9	279.8	363.6;
%     182	0	42.3	98.3	160.4	227.2	296.4	380.2;
%     233	0	52.4	110.8	174.0	241.1	310.3	394.9;
%     281	0	58.0	125.5	196.5	270.1	345.0	436.0;
%     382	0	67.6	137.7	211.2	287.4	364.8	459.1;
%     426	0	55.8	126.5	201.7	279.9	359.1	455.4;
%     521	0	74.1	148.1	225.4	304.9	385.5	484.2;
    ];

Zsrt = sortrows(Zun(1:23,:),1);           %first column is resting lengths
bpaRestL = Zsrt(:,1);
bpaNorm = bpaRestL/1000;        %Convert mm to meters also "normalizes" these lengths

pointz = [0,100,200,300,400,500,620];  %Pressure value for columns
Zfin = Zsrt(:,2:8);                    %Z values (i.e. force) are these columns

% Znorm = Zfin/(max(max(Zfin)));      %Scale force   
% P_norm = pointz/max(pointz);        %Scale pressure

%Now use curve fitting tool Fmax_fitTool.sfit
md3 = @(x,y) 0.4909*y.*atan(0.03171*y.*(x-0.0075));
z = feval(md3,bpaNorm,pointz);

%Create accessible color pallete for plotting
c6 = [0.392156862745098 0.56078431372549 1];
c5 = [0.470588235294118 0.368627450980392 0.941176470588235];
c4 = [0.862745098039216 0.149019607843137 0.498039215686275];
c3 = [0.996078431372549 0.380392156862745 0];
c2 = [1 0.690196078431373 0];
c1 = [0 0 0];

fig3d = figure;
ax1 = axes('Parent',fig3d);
hold(ax1,'on');
s1 = plot(bpaNorm(~isnan(Zfin(:,7))),Zfin(~isnan(Zfin(:,7)),7),'o',bpaNorm,z(:,7),'--');
set(s1,'LineWidth',2,'Color',c1);
s2 = plot(bpaNorm(~isnan(Zfin(:,6))),Zfin(~isnan(Zfin(:,6)),6),'o',bpaNorm,z(:,6),'--');
set(s2,'LineWidth',2,'Color',c2);
s3 = plot(bpaNorm(~isnan(Zfin(:,5))),Zfin(~isnan(Zfin(:,5)),5),'o',bpaNorm,z(:,5),'--');
set(s3,'LineWidth',2,'Color',c3);
s4 = plot(bpaNorm(~isnan(Zfin(:,4))),Zfin(~isnan(Zfin(:,4)),4),'o',bpaNorm,z(:,4),'--');
set(s4,'LineWidth',2,'Color',c4);
s5 = plot(bpaNorm(~isnan(Zfin(:,3))),Zfin(~isnan(Zfin(:,3)),3),'o',bpaNorm,z(:,3),'--');
set(s5,'LineWidth',2,'Color',c5);
s6 = plot(bpaNorm(~isnan(Zfin(:,2))),Zfin(~isnan(Zfin(:,2)),2),'o',bpaNorm,z(:,2),'--');
set(s6,'LineWidth',2,'Color',c6);
xlim(ax1,[0 1.1]);
title('F_{max} vs. l_{rest} and P_{max}, \phi10 mm')
xlabel('l_{rest}, m','FontWeight','bold')
ylabel('F_{max}, N','FontWeight','bold')
hold(ax1,'off');
set(ax1,'FontSize',12,'FontWeight','bold','LineWidth',2,'TickLength',...
    [0.02 0.05],'XMinorTick','on','YMinorTick','on');
leg = legend('Measured 620 kPa','Model 620 kPa','500 kPa',' ','400 kPa',' ','300 kPa',' ','200 kPa',' ','100 kPa',' ');  %Note, use '' w/o space to remove dash from legend, use ' ' w/ space to include dashed lines in legend
set(leg,...
    'Position',[0.753623173004846 0.257738103327298 0.223602479696274 0.565476174297787]);
% hold off

%% Fit for maximum strain
%Hint: no relationship found for max strain either as a linear or 2nd
%degree polynomial function of resting length. Hypothesis of linear
%relationship is not shown. Maximum strain may instead be dependent on
%Festo BPA batch number, age of BPA, number of full spirals in BPA length,
%being stored in kinked position, or other factors. Using a tape measure
%has +/- 1 mm accuracy which can affect the maximum strain listed.

%ft = fittype('a+b*x');
[max_k, gof2, output2] = fit(restingL',kmax','poly1','Exclude',[13],'Normalize','on')
p3 = feval(max_k,x);

figure
plot(restingL,kmax,'o',x,p3,'--')
title('\epsilon_{max} vs. Resting Length')
xlabel('Resting Length, mm')
ylabel('Maximum Strain')
legend('Data','Model Fit')