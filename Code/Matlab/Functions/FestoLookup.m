%% Create Festo force lookup tables based on their datasheets
X = linspace(0,600,7); %Pressure for interpolation

%% 40mm dia BPA
Y1 = linspace(-0.05,0.25,31);   %Contraction percent for interpolation
FestoLookup40 = zeros(size(X,2),size(Y1,2));
x40_1 = [-0.05 -0.04, -0.03, -.02 0]';
y40_1 = [2800   1900, 1000,   500 0]';   %0 kPa results
x40_2 = [-0.05 -.04  -.03  -0.01   0    0.05 0.17]';
y40_2 = [4000   3000  2000  1300  1100  500     0]';    %100 kPa results
x40_3 = [-0.05 -.03 -.01 0    0.04 0.09 0.25]';
y40_3 = [5000  3200 2500 2100 1500 1000    0]';    %200 kPa results
x40_4 = [-.05 -.03, -.02 .005 0.06 0.16, .25]';
y40_4 = [6000  4250 3900 3000 2000 1000, 300]';    %300 kPa results
x40_5 = [-.05 -.039 -.01 .005 0.05 0.11 0.20 0.25]';
y40_5 = [6000  6000 4300 4000 3000 2000 1000 600]';    %400 kPa results
x40_6 = [-.05 -.02 -.01 .008 0.04 0.09 0.15 0.23 .25]';
y40_6 = [6000 6000 5750 5000 4000 3000 2000 1000 800]';    %500 kPa results
x40_7 = [-.05 .008 0.02 0.04 .075 0.12 .175 0.25]';
y40_7 = [6000 6000 5500 5000 4000 3000 2000  990]';    %620 kPa results
FestoLookup40(1,:) = interp1(x40_1,y40_1,Y1,'pchip',0);
FestoLookup40(2,:) = interp1(x40_2,y40_2,Y1,'pchip',0);
FestoLookup40(3,:) = interp1(x40_3,y40_3,Y1,'pchip',0);
FestoLookup40(4,:) = interp1(x40_4,y40_4,Y1,'pchip',0);
FestoLookup40(5,:) = interp1(x40_5,y40_5,Y1,'pchip',0);
FestoLookup40(6,:) = interp1(x40_6,y40_6,Y1,'pchip',0);
FestoLookup40(7,:) = interp1(x40_7,y40_7,Y1,'pchip',0);

%% 20 mm BPA
Y2 = linspace(-0.04,0.25,30);   %Relative strain range for interpolation
FestoLookup20 = zeros(size(X,2),size(Y2,2));
x20_1 = [-.04 -.03 -.02 -.01 0 .25]';
y20_1 = [1000 500   260  125 0  0 ]';   %0 kPa results
x20_2 = [-.04 -.03 -.02 -.01 0   0.02 0.06 .25]';
y20_2 = [1300 800   625  400 300 125   0     0]';    %100 kPa results
x20_3 = [-.04 -.02 -.01 0.01 .05 .155 .25]';
y20_3 = [1500  875  750  450 250    0   0]';    %200 kPa results
x20_4 = [-.04 -.03 -.01 .01 .04 .12 .21 .25]';
y20_4 = [1500 1500 1000 750 500 250   0   0]';    %300 kPa results
x20_5 = [-.04 -.02    0 0.01 .04 .09 .16 .24 .25]';
y20_5 = [1500 1500 1125 1000 750 500 250   0   0]';    %400 kPa results
x20_6 = [-.04 -.01 0.01 .04 .08 .13 .195 .245 .25]';
y20_6 = [1500 1500 1250 990 750 500  250    0   0]';    %500 kPa results
x20_7 = [-.04 0.01 .035 .065 .11 .25]';
y20_7 = [1500 1500 1250 1000 750 125]';    %620 kPa results
FestoLookup20(1,:) = interp1(x20_1,y20_1,Y2,'pchip',0);
FestoLookup20(2,:) = interp1(x20_2,y20_2,Y2,'pchip',0);
FestoLookup20(3,:) = interp1(x20_3,y20_3,Y2,'pchip',0);
FestoLookup20(4,:) = interp1(x20_4,y20_4,Y2,'pchip',0);
FestoLookup20(5,:) = interp1(x20_5,y20_5,Y2,'pchip',0);
FestoLookup20(6,:) = interp1(x20_6,y20_6,Y2,'pchip',0);
FestoLookup20(7,:) = interp1(x20_7,y20_7,Y2,'pchip',0);

%% Save it
save ForceStrainTable.mat FestoLookup40 FestoLookup20

%% Plot it
figure
surf(Y1,X,FestoLookup40)
xlabel('Contraction'),ylabel('Pressure, kPA'),zlabel('Force, N')
title('40 mm BPA Force-Pressure-Contraction relationship')

figure
surf(Y2,X,FestoLookup20)
xlabel('Contraction'),ylabel('Pressure, kPA'),zlabel('Force, N')
title('20 mm BPA Force-Pressure-Contraction relationship')
