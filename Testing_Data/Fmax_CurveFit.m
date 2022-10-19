%% Look for equations of max force and max strain as functions of resting length
%  Data from experiments. See Excel spreadsheet "StraightForceTest"

%clear; clc; close all

A = [845	840	785	780	710	709	571	571	112	415	455	490	518	551	361	54	27	69	275	151	193	10	120	220	260	281	281;
    0.171597633	0.172619048	0.169426752	0.170512821	0.177464789	0.174894217	0.183887916	0.182136602	0.160714286	0.16626506	0.158241758	0.187755102	0.166023166	0.174228675	0.171745152	0.148148148	0.148148148	0.15942029	0.170909091	0.152317881	0.165803109	0.1	0.166666667	0.159090909	0.161538462	0.163701068	0.145907473;
    447.1	447.1	472	472	452.32	452.32	461.6	461.6	334	444.8222	460	451.38	455.83	453.14	436.4	238.2	135.32	271.48	412.3	347.18	396.17	8.3	343.05	377.2	377.2	418.19	402.48];


restingL = A(1,:)/1000;

kmax = A(2,:);

Fmax = A(3,:);
Fmax_norm = Fmax/max(Fmax);
%% Nonlinear fit for Fmax
modelfun = @(b,xm)b(1)*atan(b(2)*(xm-0.0075)); %based on the shape, it looks like resting length reaches a limit
%modelfun = @(b,xm)-exp(-b(1)*(xm-b(2))+b(3))+b(4);
beta0 = [301.568 20.4972];
%beta0 = [11 0.41 1.2 478];
opts = statset('fitnlm');
opts.MaxIter = 1000;
opts.Display = 'final';
opts.RobustWgtFun = 'logistic';
mdl = fitnlm(restingL,Fmax,modelfun,beta0,'Exclude',[7,27],'Options',opts)
x = linspace(0,1);
p1 = feval(mdl,x);

figure
plot(restingL,Fmax,'o',x,p1,'--')
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
bpaNorm = bpaRestL/1000;        %Convert mm to meters also normalizes these lengths

pointz = [0,100,200,300,400,500,620];  %Pressure value for columns
Zfin = Zsrt(:,2:8);                    %Z values (i.e. force) are these columns

Znorm = Zfin/(max(max(Zfin)));      %Scale force   
P_norm = pointz/max(pointz);        %Scale pressure

%Now use curve fitting tool Fmax_fitTool.sfit
 

%% Fit for maximum strain
%Hint: no relationship found for max strain either as a linear or 2nd
%degree polynomial function of resting length. Hypothesis of linear
%relationship is not shown. Maximum strain may instead be dependent on
%Festo BPA batch number, age of BPA, number of full spirals in BPA length,
%being stored in kinked position, or other factors. Using a tape measure
%has +/- 1 mm accuracy which can affect the maximum strain listed.

%ft = fittype('a+b*x');
max_k = fit(restingL',kmax','poly1','Normalize','on')
p2 = feval(max_k,x);

figure
plot(restingL,kmax,'o',x,p2,'--')
xlabel('Resting Length, mm')
ylabel('Maximum Strain')
legend('Data','Model Fit')