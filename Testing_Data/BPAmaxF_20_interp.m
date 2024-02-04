%% 20 mm Force-Pressure-Length relationship

clear;
clc;
close all;

%% 509 mm BPA resting length

RelStrain = [0,1,0.951612903000000,0.879032258000000,0.774193548000000,0.516129032000000,0.169354839000000,0.0806451610000000,0.0645161290000000,0.0564516130000000,0.0403225810000000,0.0322580650000000,0.0322580650000000,0.0241935480000000,0.0241935480000000,0.0483870970000000,0.725806452000000,0.750000000000000,0.750000000000000,0.766129032000000,0.766129032000000,0.766129032000000,0.741935484000000,0.741935484000000,0.411290323000000,0.427419355000000,0.435483871000000,0.443548387000000,0.443548387000000,0.459677419000000,0.459677419000000,0.459677419000000,0.435483871000000,0.435483871000000,0.419354839000000,0.435483871000000,0.217741935000000,0.225806452000000,0.225806452000000,0.233870968000000,0.241935484000000,0.241935484000000,0.250000000000000,0.250000000000000,0.250000000000000,0.241935484000000,0.233870968000000,0.233870968000000,0.250000000000000]';
Pressure = [0,620,500,394,295,191,95,624,600,495,403,305,202,90,240,355,335,403,500,550,600,620,447,376,187,252,300,400,501,550,601,622,525,349,201,450,130,201,299,402,450,501,603,622,553,348,252,376,476]';
Force = [0,0,0,0,0,0,0,1316,1259,1013,811,588,355,106,444,704,5,101,207,263,316,340,163,90,4,113,194,360,525,604,688,720,566,282,39,446,4.30000000000000,149,356,563,622,762,967,1006,868,456,261,515,714]';
% Force = Force*1.0208;  %Correct for a low force reading

F = scatteredInterpolant(RelStrain, Pressure, Force, 'natural','linear');
Fmax = F(0,620);
disp('Maximum force in 509 mm BPA')
disp(Fmax)

figure
title('509 mm l_{rest}')
hold on
scatter3(RelStrain, Pressure, Force)
scatter3(0,620,Fmax,25,'LineWidth',2)
hold off

E50 = RelStrain;
P50 = Pressure;
F50 = Force/Fmax;

%% 300 mm BPA resting length

RelStrain = [0,0.347222222000000,0.625000000000000,0.791666667000000,0.888888889000000,0.986111111000000,1,1,0.915824916000000,0.861952862000000,0.754208754000000,0.552188552000000,0.121212121000000,0.0133333330000000,0.0266666670000000,0.0400000000000000,0.0533333330000000,0.0800000000000000,0.0933333330000000,0.106666667000000,0.0800000000000000,0.0800000000000000,0.0666666670000000,0.386666667000000,0.400000000000000,0.400000000000000,0.413333333000000,0.413333333000000,0.400000000000000,0.400000000000000,0.386666667000000,0.386666667000000,0.746666667000000,0.746666667000000,0.746666667000000,0.746666667000000,0.746666667000000,0.746666667000000,0.746666667000000,0.746666667000000,0.746666667000000,0.986301370000000,0.986301370000000,0.972602740000000,0.972602740000000,0.958904110000000,0.958904110000000,0.986301370000000,0.986301370000000,0.986301370000000,0.602739726000000,0.602739726000000,0.616438356000000,0.616438356000000,0.616438356000000,0.616438356000000,0.630136986000000,0.616438356000000,0.602739726000000]';
Pressure = [0,206,301,399,501,605,623,622,499,395,299,200,101,136,208,304,397,502,595,620,431,354,249,297,402,499,599,621,449,350,249,200,347,449,502,545,603,625,573,474,373,621,602,571,561,552,537,581,593,615,285,250,350,405,500,600,618,550,450]';
Force = [0,0,0,0,0,0,0,0,0,0,0,0,0,223,373,575,765,970,1136,1182,783,617,390,169,342,498,656,691,423,267,108,31,35,140,193,237,291,316,267,174,76,47,38,23,19,17,9,32,38,49,55,14,135,202,312,430,450,371,258]';

F = scatteredInterpolant(RelStrain, Pressure, Force, 'natural','linear');
Fmax = F(0,620);
disp('Maximum force in 300 mm BPA')
disp(Fmax)

figure
title('300 mm l_{rest}')
hold on
scatter3(RelStrain, Pressure, Force)
scatter3(0,620,Fmax,25,'LineWidth',2)
hold off

E30 = RelStrain;
P30 = Pressure;
F30 = Force/Fmax;

%% 369 mm BPA

RelStrain = [0,1,0.0106382980000000,0.0212765960000000,0.0319148940000000,0.0425531910000000,0.0425531910000000,0.0425531910000000,0.0425531910000000,0.0319148940000000,0.0106382980000000,0.0425531910000000,0.0212765960000000,0.138297872000000,0.414893617000000,0.691489362000000,0.819148936000000,0.925531915000000,1,0.893617021000000,0.893617021000000,0.882978723000000,0.308510638000000,0.287234043000000,0.265957447000000,0.255319149000000,0.276595745000000,0.319148936000000,0.340425532000000,0.329787234000000,0.308510638000000,0.159574468000000,0.159574468000000,0.202127660000000,0.202127660000000,0.202127660000000,0.191489362000000,0.170212766000000,0.191489362000000,0.510638298000000,0.510638298000000,0.521276596000000,0.521276596000000,0.521276596000000,0.521276596000000,0.531914894000000]';
Pressure = [0,620,106,199,295,395,497,598,624,279,96,391,198,110,198,298,397,498,625,400,495,618,626,397,202,157,299,499,625,601,391,105,302,493,625,592,381,190,391,216,399,623,598,297,203,498]';
Force = [0,0,150,354,570,790,1011,1232,1289,526,127,775,348,0,0,0,0,0,0,12.5000000000000,81,163,844,455,113,26,277,627,840,796,442,13.5000000000000,414,788,1043,980,567,185,591,13,271,573,540,139,15,410]';

F = scatteredInterpolant(RelStrain, Pressure, Force, 'natural','linear');
Fmax = F(0,620);
disp('Maximum force in 369 mm BPA')
disp(Fmax)

figure
title('369 mm l_{rest}')
hold on
scatter3(RelStrain, Pressure, Force)
scatter3(0,620,Fmax,25,'LineWidth',2)
hold off

E37 = RelStrain;
P37 = Pressure;
F37 = Force/Fmax;

%% 451 mm BPA

RelStrain = [0,1,0.948717949000000,0.905982906000000,0.777777778000000,0.658119658000000,0.256410256000000,0.0940170940000000,0.0940170940000000,0.0854700850000000,0.0769230770000000,0.0683760680000000,0.0598290600000000,0.0683760680000000,0.0769230770000000,0.0769230770000000,0.0940170940000000,0.0769230770000000,0.0341880340000000,0.0427350430000000,0.0512820510000000,0.0512820510000000,0.0598290600000000,0.0598290600000000,0.0683760680000000,0.0598290600000000,0.0512820510000000,0.0341880340000000,0.119658120000000,0.162393162000000,0.170940171000000,0.170940171000000,0.179487179000000,0.188034188000000,0.188034188000000,0.188034188000000,0.179487179000000,0.170940171000000,0.162393162000000,0.316239316000000,0.324786325000000,0.341880342000000,0.341880342000000,0.350427350000000,0.350427350000000,0.358974359000000,0.341880342000000,0.324786325000000,0.316239316000000,0.705357143000000,0.696428571000000,0.687500000000000,0.687500000000000,0.687500000000000,0.696428571000000,0.705357143000000,0.705357143000000,0.982142857000000,0.973214286000000,0.973214286000000,0.973214286000000,0.973214286000000]';
Pressure = [0,620,509,398,297,205,100,623,595,501,403,198,115,251,350,448,543,290,234,360,449,502,552,604,620,469,285,164,175,249,301,402,498,598,615,555,455,346,201,250,350,451,501,550,603,621,399,300,199,612,500,400,301,348,450,550,628,615,594,571,582,600]';
Force1 = [0,0,0,0,0,0,0,1307.64480000000,1246.39680000000,1026.92480000000,802.348800000000,331.760000000000,140.870400000000,453.235200000000,676.790400000000,900.345600000000,1119.81760000000,538.982400000000,472.630400000000,764.579200000000,971.801600000000,1094.29760000000,1209.64800000000,1319.89440000000,1351.53920000000,1001.40480000000,572.668800000000,290.928000000000,110.246400000000,272.553600000000,387.904000000000,603.292800000000,810.515200000000,1023.86240000000,1059.59040000000,934.032000000000,724.768000000000,491.004800000000,183.744000000000,136.787200000000,323.593600000000,510.400000000000,597.168000000000,683.936000000000,773.766400000000,806.432000000000,425.673600000000,255.200000000000,76.5600000000000]';
% Force1 = Force1*1.0208; %Correct for a low force reading
Force2 = [330,221,114,11,62,172,278,361,27,18,9,19,31]';
Force = [Force1; Force2];

F = scatteredInterpolant(RelStrain, Pressure, Force, 'natural','linear');
Fmax = F(0,620);
disp('Maximum force in 451 mm BPA')
disp(Fmax)

figure
title('451 mm l_{rest}')
hold on
scatter3(RelStrain, Pressure, Force)
scatter3(0,620,Fmax,25,'LineWidth',2)
hold off

E45 = RelStrain;
P45 = Pressure;
F45 = Force/Fmax;

%% Plot normalized values together, including from Festo
load FestoData.mat XX20 YY20 ZZ20

YY20 = YY20*620;

figure
hold on
scatter3(E50,P50,F50,'DisplayName','509 mm l_{rest}')
scatter3(E45,P45,F45,'DisplayName','451 mm l_{rest}')
scatter3(E37,P37,F37,'DisplayName','369 mm l_{rest}')
scatter3(E30,P30,F30,'DisplayName','300 mm l_{rest}')
scatter3(XX20,YY20,ZZ20,[],'d','DisplayName','Festo')
hold off
lgd = legend;

%% Combine experimental data for fitting tool

Ax = [E50; E45; E37; E30];
Ay = [P50; P45; P37; P30]/620;
Az = [F50; F45; F37; F30];

save data20mm.mat Ax Ay Az