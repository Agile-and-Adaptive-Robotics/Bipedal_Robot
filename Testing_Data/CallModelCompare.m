%% Call the Model Comparison 2 function


clear; clc; close all

load Data4compare10.mat W1 W2
load datacompare20.mat E30 P30 F30 E45 P45 F45

D = {10; 20}; %Diameters

Lrest = {
    [0.257 0.233] %10 diameter top row
    [0.300 0.450] %20 mm diameter bottom row
};

Lmin = {
    [0.221 0.193];
    [0.225 0.335]
};

M30 = [E30, P30/620, F30]; M45 = [E45, P45/620, F45];

Data = {
    {W1, W2}
    {M30, M45}
};

set(groot, ...
    'defaultAxesLineWidth',3, ...
    'defaultAxesFontSize',12, ...
    'defaultAxesFontWeight','bold')
Goof = ModelComparison2(D, Lrest, Lmin, Data, 'MartensBest','on','SarosiBest','on')