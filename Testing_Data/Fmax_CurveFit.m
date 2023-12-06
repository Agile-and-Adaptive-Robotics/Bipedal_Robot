%% Look for equations of max force and max strain as functions of resting length
%  Data from experiments. See Excel spreadsheet "StraightForceTest"

clear; clc; close all

A = [845	840	785	780	710	709	571	571	112	415	455	490	518	551	361	54	27	69	275	151	193	10	120	220	260	281	281;
    0.171597633	0.172619048	0.169426752	0.170512821	0.177464789	0.174894217	0.183887916	0.182136602	0.160714286	0.16626506	0.158241758	0.187755102	0.166023166	0.174228675	0.171745152	0.148148148	0.148148148	0.15942029	0.170909091	0.152317881	0.165803109	0.1	0.166666667	0.159090909	0.161538462	0.163701068	0.145907473;
    447.1	447.1	472	472	452.32	452.32	461.6	461.6	334	444.8222	460	451.38	455.83	453.14	436.4	238.2	135.32	271.48	412.3	347.18	396.17	8.3	343.05	377.2	377.2	418.19	402.48];


restingL = A(1,:)/1000;  %Resting length

kmax = A(2,:);          %Maximum strain

Fmax = A(3,:);          %Maximum force
Fmax_norm = Fmax/max(Fmax);
%% Nonlinear fit for Fmax
modelfun = @(b,xm)b(1)*620*atan(b(2)*(xm-0.0075)*620); %based on the shape, it looks like resting length reaches a limit
%modelfun = @(b,xm)-exp(-b(1)*(xm-b(2))+b(3))+b(4);
beta0 = [0.4864 0.03306];
%beta0 = [11 0.41 1.2 478];
opts = statset('fitnlm');
opts.MaxIter = 1000;
opts.Display = 'final';
opts.RobustWgtFun = 'logistic';
mdl = fitnlm(restingL,Fmax,modelfun,beta0,'Exclude',[7,27],'Options',opts)
x = linspace(0,1);
p1 = feval(mdl,x);

md2 = @(x) 0.4864*620*atan(0.03306*620*(x-0.0075));     %from Fmax_fitTool
p2 = feval(md2,x);

figure
plot(restingL,Fmax,'o',x,p2,'--')
title('Max. Force vs. Resting Length')
xlabel('Resting Length, mm')
ylabel('Maximum Force, N')
legend('Data','Model Fit')

%% Fit maximum force (i.e. zero contraction BPA) for different resting lengths and pressures
%Load data from StraightForceTest Fmax(L) tab. It uses data from some of
%Lawrence's tests.

% Force (z-axis) unsorted
Zun = [840	0	NaN	114.32	191.8	275.2	355.4	447.1
        780	0	NaN	139.2	216	299.7	378.1	472
        709	0	NaN	124.21	201.4	276.02	355.9	452.32
        571	0	NaN	136.1	212.39	288.2	366.94	461.6
        112	0	NaN	NaN	NaN	NaN	NaN	334
        415	0	NaN	NaN	NaN	NaN	NaN	444.8
        455	0	NaN	NaN	NaN	NaN	NaN	460
        490	0	NaN	NaN	NaN	NaN	NaN	451.38
        518	0	NaN	NaN	NaN	NaN	NaN	455.83
        551	0	NaN	129.26	201.3	282.34	360.24	453.14
        361	0	NaN	121.5	192.91	269	344.5	436.4
        54	0	NaN	56.1	103.37	155.45	207	266.5
        27	0	NaN	NaN	NaN	NaN	NaN	135.32
        69	0	NaN	NaN	NaN	NaN	NaN	271.48
        275	0	NaN	NaN	NaN	NaN	NaN	412.3
        151	0	NaN	NaN	NaN	NaN	NaN	347.11
        193	0	NaN	NaN	NaN	NaN	NaN	396.17
        10	0	NaN	NaN	NaN	NaN	NaN	8.3
        120	0	19.7	66	124.82	195.1	258.3	341.5
        220	0	23	80.8	146.2	215.6	289.8	377.23
        260	0	31.9	89.7	154.6	225.3	294.8	385.4
        281	0	54.65	117.8	185	257.2	329.9	415.8
        281	0	37.1	101.1	170.4	242.3	314.2	406];

Zsrt = sortrows(Zun);           %first column is resting lengths
bpaRestL = Zsrt(:,1);
bpaNorm = bpaRestL/1000;        %Convert mm to meters also "normalizes" these lengths

pointz = [0,100,200,300,400,500,620];  %Pressure value for columns
Zfin = Zsrt(:,2:8);                    %Z values (i.e. force) are these columns

% Znorm = Zfin/(max(max(Zfin)));      %Scale force   
% P_norm = pointz/max(pointz);        %Scale pressure

%Now use curve fitting tool Fmax_fitTool.sfit
md3 = @(x,y) 0.4864*y.*atan(0.03306*y.*(x-0.0075));
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
[max_k, gof2, output2] = fit(restingL',kmax','poly1','Exclude',[22],'Normalize','on')
p3 = feval(max_k,x);

figure
plot(restingL,kmax,'o',x,p3,'--')
title('Max. Strain vs. Resting Length')
xlabel('Resting Length, mm')
ylabel('Maximum Strain')
legend('Data','Model Fit')