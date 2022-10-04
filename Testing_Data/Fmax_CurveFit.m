%% Look for equations of max force and max strain as functions of resting length
%  Data from experiments. See Excel spreadsheet "StraightForceTest"

clear; clc; close all

restingL = [845	840	785	780	710	709	571	571	112	415	455	490	518	551	361	54	27	69	275	151	193 10];

kmax = [0.171597633	0.172619048	0.169426752	0.170512821	0.177464789	0.174894217	0.183887916	0.182136602	0.160714286	0.16626506	0.158241758	0.187755102	0.166023166	0.174228675	0.171745152	0.148148148	0.148148148	0.15942029	0.170909091	0.152317881	0.165803109 0.1];

Fmax = [447.1	447.1	472	472	452.32	452.32	461.6	461.6	334	445	464.5	458.6	490	453.14	436.4	238.2	135.32	271.48	412.3	347.18	396.17  8.3];

%% Nonlinear fit for Fmax
modelfun = @(b,x)b(1)*tanh(b(2)*x+b(3)); %based on the shape, it looks like resting length reaches a limit
beta0 = [450.9 0.0089705 0];
opts = statset('fitnlm');
opts.MaxIter = 1000;
mdl = fitnlm(restingL,Fmax,modelfun,beta0,'Exclude',[1,3,5,7],'Options',opts)
x = linspace(0,1000);
p1 = feval(mdl,x);

figure
plot(restingL,Fmax,'o',x,p1,'--')
xlabel('Resting Length, mm')
ylabel('Maximum Force, N')
legend('Data','Model Fit')

%% Fit for maximum strain

%ft = fittype('a+b*x');
max_k = fit(restingL',kmax','poly2','Normalize','on')
p2 = feval(max_k,x);

figure
plot(restingL,kmax,'o',x,p2,'--')
xlabel('Resting Length, mm')
ylabel('Maximum Strain')
legend('Data','Model Fit')