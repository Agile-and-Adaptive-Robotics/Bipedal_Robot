clear;
clc;
close all;

%% Create Festo force lookup tables based on their datasheets
Y = linspace(0,600,7);                      %Pressure for interpolation

%% 40mm dia BPA
X1 = linspace(-0.05,0.25,31);   %Contraction percent for interpolation
FestoLookup40 = zeros(size(Y,2),size(X1,2));
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


x40 = [x40_1; x40_2; x40_3; x40_4; x40_5; x40_6; x40_7];
y40 = [0*ones(length(x40_1),1); 
       100*ones(length(x40_2),1); 
       200*ones(length(x40_3),1);
       300*ones(length(x40_4),1);
       400*ones(length(x40_5),1);
       500*ones(length(x40_6),1);
       600*ones(length(x40_7),1)];
z40 = [z40_1;
       z40_2;
       z40_3;
       z40_4;
       z40_5;
       z40_6;
       z40_7];


%% 20 mm BPA
X2 = linspace(-0.04,0.25,30);   %Strain range for interpolation
FestoLookup20 = zeros(size(Y,2),size(X2,2));
x20_1 = [-.04 -.03 -.02 -.01 0 0.03 .125 .25]';
z20_1 = [ 684  416  229   96 0 -160 -270 -287]';                   %0 kPa force, N
x20_2 = [-.04 -.03 -.02 -.01 0   0.02 0.06 0.17 .25]';
z20_2 = [ 940  672  486  352 253 123    -8 -147 -218]';       %100 kPa force, N
x20_3 = [-.04 -.02 -.01 0.01 .05 .12 .155 .25]';
z20_3 = [1194  743  608  429 240 112 8.7  -150]';            %200 kPa force, N
x20_4 = [-.04 -.03 -.01 .01 .04 .12 .21 .25]';
z20_4 = [1451 1187  864 678 507 241  13 -84]';          %300 kPa force, N
x20_5 = [-.04 -.02    0 .01 .04 .09 .16 .24 .25]';
z20_5 = [1707 1257 1012 926 737 518 271  13 -19]';     %400 kPa force, N
x20_6 = [-.04 -.02    0 0.01 .04 .08 .13 .195 .22 .25]';
z20_6 = [1963 1514 1264 1173 966 756 531  262 162  45]';    %500 kPa force, N
x20_7 = [-.04 0.00 .035 .065 .11 .25]';
z20_7 = [2219 1515 1228 1040 791 107]';                 %600 kPa force, N

x20 = [x20_1; x20_2; x20_3; x20_4; x20_5; x20_6; x20_7];
y20 = [0*ones(length(x20_1),1); 
       100*ones(length(x20_2),1); 
       200*ones(length(x20_3),1);
       300*ones(length(x20_4),1);
       400*ones(length(x20_5),1);
       500*ones(length(x20_6),1);
       600*ones(length(x20_7),1)];
z20 = [z20_1; z20_2; z20_3; z20_4; z20_5; z20_6; z20_7];

%% From Curve Fit app generate code

% Initialize arrays to store fits and goodness-of-fit.
fitresult = cell( 4, 1 );
gof = struct( 'sse', cell( 7, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit: 'Exponential 40'.
[xData, yData, zData] = prepareSurfaceData( x40, y40, z40 );

% Set up fittype and options.
ft = fittype( 'a0+exp(-a1.*x+a2)+y.*exp(-a3.*x+a4)+a5*y', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.655740699156587 0.0357116785741896 0.849129305868777 1.2 0.6787 0.757740130578333];

% Fit model to data.
[fitresult{1}, gof(1)] = fit( [xData, yData], zData, ft, opts );
f40 = fitresult{1}
gof40 = gof(1)

FestoLookup40(1,:) = f40(X1,0);
FestoLookup40(2,:) = f40(X1,100);
FestoLookup40(3,:) = f40(X1,200);
FestoLookup40(4,:) = f40(X1,300);
FestoLookup40(5,:) = f40(X1,400);
FestoLookup40(6,:) = f40(X1,500);
FestoLookup40(7,:) = f40(X1,600);

% Create a figure for the plots.
figure( 'Name', 'Exponential 40' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult{1}, [xData, yData], zData );
legend( h, 'Exponential 40', 'z40 vs. x40, y40', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'x40', 'Interpreter', 'latex' );
ylabel( 'y40', 'Interpreter', 'latex' );
zlabel( 'z40', 'Interpreter', 'latex' );
grid on
view( -0.9, 1.7 );

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult{1}, [xData, yData], zData, 'Style', 'Residual' );
legend( h, 'Exponential 40 - residuals', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'x40', 'Interpreter', 'latex' );
ylabel( 'y40', 'Interpreter', 'latex' );
zlabel( 'z40', 'Interpreter', 'latex' );
grid on
view( -0.9, 1.7 );

%% Fit: 'Exponential 20'.
[xData, yData, zData] = prepareSurfaceData( x20, y20, z20 );

% Set up fittype and options.
ft = fittype( 'a0+exp(-a1*x+a2)+y.*exp(-(a3*x).^2+a4)+a5*y', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxFunEvals = 6000;
opts.MaxIter = 4000;
opts.StartPoint = [0.530797553008973 0.779167230102011 0.934010684229183 0.12990620847373 0.5688 -0.4694];

% Fit model to data.
[fitresult{2}, gof(2)] = fit( [xData, yData], zData, ft, opts );
f20 = fitresult{2}
gof20 = gof(2)

FestoLookup20(1,:) = f20(X2,0);
FestoLookup20(2,:) = f20(X2,100);
FestoLookup20(3,:) = f20(X2,200);
FestoLookup20(4,:) = f20(X2,300);
FestoLookup20(5,:) = f20(X2,400);
FestoLookup20(6,:) = f20(X2,500);
FestoLookup20(7,:) = f20(X2,600);

% Create a figure for the plots.
figure( 'Name', 'Exponential 20' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult{2}, [xData, yData], zData);
legend( h, 'Exponential 20', 'z20 vs. x20, y20', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'x20', 'Interpreter', 'latex' );
ylabel( 'y20', 'Interpreter', 'latex' );
zlabel( 'z20', 'Interpreter', 'latex' );
grid off
view( -12.3, 35.2 );

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult{2}, [xData, yData], zData, 'Style', 'Residual' );
legend( h, 'Exponential 20 - residuals', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'x20', 'Interpreter', 'latex' );
ylabel( 'y20', 'Interpreter', 'latex' );
zlabel( 'z20', 'Interpreter', 'latex' );
grid off
view( -12.3, 35.2 );

%% Fit: 'Polynomial 20'.
[xData, yData, zData] = prepareSurfaceData( x20, y20, z20 );

% Set up fittype and options.
ft = fittype( 'poly51' );

% Fit model to data.
[fitresult{3}, gof(3)] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Create a figure for the plots.
figure( 'Name', 'Polynomial 20' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult{3}, [xData, yData], zData );
legend( h, 'Polynomial 20', 'z20 vs. x20, y20', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'x20', 'Interpreter', 'latex' );
ylabel( 'y20', 'Interpreter', 'latex' );
zlabel( 'z20', 'Interpreter', 'latex' );
grid on
view( -29.8, 0.9 );

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult{3}, [xData, yData], zData, 'Style', 'Residual' );
legend( h, 'Polynomial 20 - residuals', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'x20', 'Interpreter', 'latex' );
ylabel( 'y20', 'Interpreter', 'latex' );
zlabel( 'z20', 'Interpreter', 'latex' );
grid on
view( -29.8, 0.9 );

%% Fit: 'Polynomial 40'.
[xData, yData, zData] = prepareSurfaceData( x40, y40, z40 );

% Set up fittype and options.
ft = fittype( 'poly51' );

% Fit model to data.
[fitresult{4}, gof(4)] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

% Create a figure for the plots.
figure( 'Name', 'Polynomial 40' );

% Plot fit with data.
subplot( 2, 1, 1 );
h = plot( fitresult{4}, [xData, yData], zData );
legend( h, 'Polynomial 40', 'z40 vs. x40, y40', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'x40', 'Interpreter', 'latex' );
ylabel( 'y40', 'Interpreter', 'latex' );
zlabel( 'z40', 'Interpreter', 'latex' );
grid off
view( 1.4, 7.2 );

% Plot residuals.
subplot( 2, 1, 2 );
h = plot( fitresult{4}, [xData, yData], zData, 'Style', 'Residual' );
legend( h, 'Polynomial 40 - residuals', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'x40', 'Interpreter', 'latex' );
ylabel( 'y40', 'Interpreter', 'latex' );
zlabel( 'z40', 'Interpreter', 'latex' );
grid off
view( 1.4, 7.2 );


%% Save it
save FestoLookup.mat f40 f20

%% Plot it
figure
plot(f40, [x40, y40],z40)
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Pressure, $kPA$','interpreter','latex'),zlabel('\bf Force, $N$','interpreter','latex')
title('\bf 40 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')

figure
plot(f20, [x20, y20],z20)
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Pressure, $kPA$','interpreter','latex'),zlabel('\bf Force, $N$','interpreter','latex')
title('\bf 40 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')

figure
surf(X1,Y,FestoLookup40)
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Pressure, $kPA$','interpreter','latex'),zlabel('\bf Force, $N$','interpreter','latex')
title('40 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')

figure
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Force, $N$','interpreter','latex')
title('\bf 40 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')
hold on
plot(X1,FestoLookup40(7,:),'DisplayName','600 kPa')
plot(X1,FestoLookup40(6,:),'DisplayName','500 kPa')
plot(X1,FestoLookup40(5,:),'DisplayName','400 kPa')
plot(X1,FestoLookup40(4,:),'DisplayName','300 kPa')
plot(X1,FestoLookup40(3,:),'DisplayName','200 kPa')
plot(X1,FestoLookup40(2,:),'DisplayName','100 kPa')
plot(X1,FestoLookup40(1,:),'DisplayName','  0 kPa')
lgd40 = legend('interpreter','latex');
title(lgd40,'\bf Pressure')
hold off

figure
surf(X2,Y,FestoLookup20)
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Pressure, kPA','interpreter','latex'),zlabel('\bf Force, N','interpreter','latex')
title('\bf 20 $ mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')

figure
xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Force, $N$','interpreter','latex')
title('\bf 20 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')
hold on
plot(X2,FestoLookup20(7,:),'DisplayName','600 kPa')
plot(X2,FestoLookup20(6,:),'DisplayName','500 kPa')
plot(X2,FestoLookup20(5,:),'DisplayName','400 kPa')
plot(X2,FestoLookup20(4,:),'DisplayName','300 kPa')
plot(X2,FestoLookup20(3,:),'DisplayName','200 kPa')
plot(X2,FestoLookup20(2,:),'DisplayName','100 kPa')
plot(X2,FestoLookup20(1,:),'DisplayName','  0 kPa')
lgd40 = legend('interpreter','latex');
title(lgd40,'\bf Pressure')
hold off
