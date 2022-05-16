% clear;
% clc;
% close all;

%% Create Festo force lookup tables based on their datasheets
X = linspace(0,600,7);                      %Pressure for interpolation

%% 40mm dia BPA
Y1 = linspace(-0.05,0.25,31);   %Contraction percent for interpolation
FestoLookup40 = zeros(size(X,2),size(Y1,2));
x40_1 = [-0.05 -0.04 -0.03, -.02   -.01  0 0.25]';
z40_1 = [2141   1325 789.6   427.4 176.6 0 -383]';                  %0 kPa force, N
x40_2 = [-0.05  -.04  -.03  -0.02 -0.01     0   0.05 0.12 0.18 0.25]';
z40_2 = [3317   2482  1922   1532  1250  1040   500  193    22 -155]';    %100 kPa force, N
x40_3 = [-0.05 -.03 -.01    0 0.02  0.04 0.08 0.12 0.20 0.25]';
z40_3 = [4489  3053 2321 2078 1720  1463 1093  809  340   76]';         %200 kPa force, N
x40_4 = [-.05 -.03 -.02 .005 0.03 0.06 0.16, .25]';
z40_4 = [5657 4179 3733 3000 2513 2084 1000, 311]';         %300 kPa force, N
x40_5 = [-0.05 -.041 -.01 .005 0.05 0.11 0.20 0.25]';
z40_5 = [ 6817  6000 4454 4009 3074 2178 1091 551]';     %400 kPa force, N
x40_6 = [-.05 -.022 -.01 .008 0.04 0.09 0.15 0.23 .25]';
z40_6 = [7967  6006 5513 4933 4142 3173 2200 1062 800]';    %500 kPa force, N
x40_7 = [-.05 -.02 .006 0.02 0.04 .075 0.12 .175 0.25]';
z40_7 = [9105 7000 6000 5568 5031 4214 3293 2286 1051]';     %600 kPa force, N
FestoLookup40(1,:) = interp1(x40_1,z40_1,Y1,'linear');
FestoLookup40(2,:) = interp1(x40_2,z40_2,Y1,'linear');
FestoLookup40(3,:) = interp1(x40_3,z40_3,Y1,'linear');
FestoLookup40(4,:) = interp1(x40_4,z40_4,Y1,'linear');
FestoLookup40(5,:) = interp1(x40_5,z40_5,Y1,'linear');
FestoLookup40(6,:) = interp1(x40_6,z40_6,Y1,'linear');
FestoLookup40(7,:) = interp1(x40_7,z40_7,Y1,'linear');

x40 = [x40_1; x40_2; x40_3; x40_4; x40_5; x40_6; x40_7];
y40 = [0*ones(length(x40_1),1); 
       100*ones(length(x40_2),1); 
       200*ones(length(x40_3),1);
       300*ones(length(x40_4),1);
       400*ones(length(x40_5),1);
       500*ones(length(x40_6),1);
       600*ones(length(x40_7),1)];
z40 = [z40_1; z40_2; z40_3; z40_4; z40_5; z40_6; z40_7];
[f40, gof40] = fit([x40, y40],z40,'poly51','Normalize','on')

figure
hold on
plot( f40, [x40, y40],z40)

hold off

%% 20 mm BPA
Y2 = linspace(-0.04,0.25,30);   %Relative strain range for interpolation
FestoLookup20 = zeros(size(X,2),size(Y2,2));
x20_1 = [-.04 -.03 -.02 -.01 0 .25]';
z20_1 = [ 684 500   260  125 0  0 ]';                   %0 kPa force, N
x20_2 = [-.04 -.03 -.02 -.01 0   0.02 0.06 .25]';
z20_2 = [1300 800   625  400 300 125   0     0]';       %100 kPa force, N
x20_3 = [-.04 -.02 -.01 0.01 .05 .155 .25]';
z20_3 = [1500  875  750  450 250    0   0]';            %200 kPa force, N
x20_4 = [-.04 -.03 -.01 .01 .04 .12 .21 .25]';
z20_4 = [1500 1500 1000 750 500 250   0   0]';          %300 kPa force, N
x20_5 = [-.04 -.02    0 0.01 .04 .09 .16 .24 .25]';
z20_5 = [1500 1500 1125 1000 750 500 250   0   0]';     %400 kPa force, N
x20_6 = [-.04 -.01 0.01 .04 .08 .13 .195 .245 .25]';
z20_6 = [1500 1500 1250 990 750 500  250    0   0]';    %500 kPa force, N
x20_7 = [-.04 0.01 .035 .065 .11 .25]';
z20_7 = [1500 1500 1250 1000 750 125]';                 %620 kPa force, N
FestoLookup20(1,:) = interp1(x20_1,z20_1,Y2,'linear');
FestoLookup20(2,:) = interp1(x20_2,z20_2,Y2,'linear');
FestoLookup20(3,:) = interp1(x20_3,z20_3,Y2,'linear');
FestoLookup20(4,:) = interp1(x20_4,z20_4,Y2,'linear');
FestoLookup20(5,:) = interp1(x20_5,z20_5,Y2,'linear');
FestoLookup20(6,:) = interp1(x20_6,z20_6,Y2,'linear');
FestoLookup20(7,:) = interp1(x20_7,z20_7,Y2,'linear');

%% Save it
save FestoLookup.mat FestoLookup40 FestoLookup20 f40 f20

%% Plot it
figure
surf(Y1,X,FestoLookup40)
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Pressure, $kPA$','interpreter','latex'),zlabel('\bf Force, $N$','interpreter','latex')
title('40 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')

figure
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Force, $N$','interpreter','latex')
title('\bf 40 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')
hold on
plot(Y1,FestoLookup40(7,:),'DisplayName','600 kPa')
plot(Y1,FestoLookup40(6,:),'DisplayName','500 kPa')
plot(Y1,FestoLookup40(5,:),'DisplayName','400 kPa')
plot(Y1,FestoLookup40(4,:),'DisplayName','300 kPa')
plot(Y1,FestoLookup40(3,:),'DisplayName','200 kPa')
plot(Y1,FestoLookup40(2,:),'DisplayName','100 kPa')
plot(Y1,FestoLookup40(1,:),'DisplayName','  0 kPa')
lgd40 = legend('interpreter','latex');
title(lgd40,'\bf Pressure')
hold off

figure
surf(Y2,X,FestoLookup20)
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Pressure, kPA','interpreter','latex'),zlabel('\bf Force, N','interpreter','latex')
title('\bf 20 $ mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')
