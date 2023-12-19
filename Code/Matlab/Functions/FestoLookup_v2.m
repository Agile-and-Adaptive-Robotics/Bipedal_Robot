clear;
clc;
close all;

%% Create Festo force lookup tables based on their datasheets
Y1 = [0, 100, 200, 300, 400, 500, 600, 620];                %Pressure for interpolation, 20 & 40 mm
Y2 = [0, 100, 200, 300, 400, 500, 600, 620, 700, 800];      %Pressure for interpolation, 10 mm

z10max = 488.4;         %Max force in 10mm BPA @ 620 kPa
z10maxFesto = 630;      %Max force in 10mm BPA, specified by Festo

z20max = 1565;          %Max force in 20mm BPA @ 620 kPa
z20maxFesto = 1500;     %Max force in 20mm BPA, specified by Festo

z40max = 6398.4;        %Max force in 40mm BPA @ 620 kPa
z40maxFesto = 6000;     %Max force in 40mm BPA, specified by Festo

Pmax = 620;             %Max pressure in 20 & 40 mm BPAs (festo limit is 600 kPa)
P10max = 620;           %Max pressure in 10 mm BPAs (festo limit is 800 kPa)

E10max = 0.2238;        %Max contraction in 10 mm BPAs @ 620 kPa
E20max = 0.275;         %Max contraction in 20 mm BPA @ 620 kPa
E40max = 0.32;          %Max contraction in 40 mm BPAs at 620 kPa


% %disregard sretching more than Festo limit
% Emin10= -0.03;
% Emin20= -0.04;
% Emin40= -0.05;

%disregard sretching more than 2% resting length
Emin10= -0.02;
Emin20= -0.02;
Emin40= -0.02;

X1 = linspace(-0.02,0.32,35);   %Strain range for interpolation, 40mm
X2 = linspace(-0.02,0.25,28);   %Strain range for interpolation, 10 & 20mm
%% 40mm dia BPA
x40 = cell(length(Y1),1);
z40 = cell(length(x40),1);

x40{1} = [-0.05 -0.04 -0.03, -.02   -.01  0 0.25]';
z40{1} = [2141   1325 789.6   427.4 176.6 0 -383]';                  %0 kPa force, N
x40{2} = [-0.05  -.04  -.03  -0.02 -0.01     0   0.05 0.12 0.18 0.25]';
z40{2} = [3317   2482  1922   1532  1250  1040   500  193    22 -155]';    %100 kPa force, N
x40{3} = [-0.05 -.03 -.01    0 0.02  0.04 0.08 0.12 0.20 0.25 0.26]';
z40{3} = [4489  3053 2321 2078 1720  1463 1093  809  340   76 24.7]';         %200 kPa force, N
x40{4} = [-.05 -.03 -.02 .005 0.03 0.06 0.16, .25 0.29]';
z40{4} = [5657 4179 3733 3000 2513 2084 1000, 311 4.2]';         %300 kPa force, N
x40{5} = [-0.05 -.041 -.01 .005 0.05 0.11 0.20 0.25 0.3]';
z40{5} = [ 6817  6000 4454 4009 3074 2178 1091 551 41.5]';     %400 kPa force, N
x40{6} = [-.05  -.02 -.01 .008 0.04 0.09 0.15 0.23 .25 0.31]';
z40{6} = [7967  5915 5513 4933 4142 3173 2200 1062 800 37.6]';    %500 kPa force, N
x40{7} = [-.05 -.02 .006 0.02 0.04 .075 0.12 .175 0.25 0.31]';
z40{7} = [9105 7000 6000 5568 5031 4214 3293 2286 1051 141.2]';     %600 kPa force, N
x40{8} = [-.05 -.02 .013 0.02 0.04 .075 0.12 .175 0.25 0.32]';
z40{8} = [9331 7210 5970 5759 5209 4368 3418 2378 1102 12]';     %620 kPa force, N

% x40 = [x40{1}; x40{2}; x40{3}; x40{4}; x40{5}; x40{6}; x40{7}; x40{8}];
y40 = cell(length(x40),1);
    for i = 1:length(x40)
        y40{i} = [Y1(i).*ones(length(x40{i}),1)]; 
    end

%Remove values with a lot of stretching or those above the Festo specified
%value
Festemp = cell2mat([x40, y40, z40]);
Fest40 = Festemp((Festemp(:,1)>=Emin40 &Festemp(:,3)<=z40maxFesto &Festemp(:,3)>=0),:);

x_40 = Fest40(:,1); y_40=Fest40(:,2); z_40=Fest40(:,3);     %convert to column vectors again


x40norm = x_40/E40max;
y40norm = y_40/Pmax;
z40norm = z_40/z40max;

[XX40, YY40, ZZ40] = prepareSurfaceData(x40norm, y40norm, z40norm);
%% 20 mm BPA
x20 = cell(length(Y1),1);
z20 = cell(length(x20),1);

x20{1} = [-.04 -.03 -.02 -.01 0 0.03 .125 .25]';
z20{1} = [ 684  416  229   96 0 -160 -270 -287]';                   %0 kPa force, N
x20{2} = [-.04 -.03 -.02 -.01 0   0.02 0.06 0.17 .25]';
z20{2} = [ 940  672  486  352 253 123    -8 -147 -218]';       %100 kPa force, N
x20{3} = [-.04 -.02 -.01 0.01 .05 .12 .155 .25]';
z20{3} = [1194  743  608  429 240 112 8.7  -150]';            %200 kPa force, N
x20{4} = [-.04 -.03 -.01 .01 .04 .12 .215 .25]';
z20{4} = [1451 1187  864 678 507 241  0.5 -84]';          %300 kPa force, N
x20{5} = [-.04 -.02    0 .01 .04 .09 .16 .24 .25]';
z20{5} = [1707 1257 1012 926 737 518 271  13 -19]';     %400 kPa force, N
x20{6} = [-.04 -.019    0 0.01 .04 .08 .13 .195 .22 .25  .26]';
z20{6} = [1963  1498 1264 1173 966 756 531  262 162  45  5.8]';    %500 kPa force, N
x20{7} = [-.04 0.00 .002 .035 .065 .11 .25]';
z20{7} = [2219 1515 1495 1228 1040 791 107]';                 %600 kPa force, N
x20{8} = [-.04 0.00 .007 .035 .065 .11 .25]';
z20{8} = [2270 1565 1497 1274 1082 826 119]';                 %620 kPa force, N


y20 = cell(length(x20),1);
    for i = 1:length(x20)
        y20{i} = [Y1(i).*ones(length(x20{i}),1)]; 
    end

%Remove values with a lot of stretching or those above the Festo specified
%value
Festemp = cell2mat([x20, y20, z20]);
Fest20 = Festemp((Festemp(:,1)>=Emin20 &Festemp(:,3)<=z20maxFesto &Festemp(:,3)>=0),:);

x_20 = Fest20(:,1); y_20=Fest20(:,2); z_20=Fest20(:,3);     %convert to column vectors again

x20norm = x_20/E20max;
y20norm = y_20/Pmax;
z20norm = z_20/z20max;

[XX20, YY20, ZZ20] = prepareSurfaceData(x20norm, y20norm, z20norm);
%% 10 mm BPA
X3 = linspace(-0.03,0.25,29);   %Strain range for interpolation
FestoLookup10 = zeros(size(Y2,2),size(X3,2));

x10 = cell(length(Y2),1);
z10 = cell(length(x10),1);

x10{1} = [-.03   -.025  -.02   -.015  -.01  -.005 0]';
z10{1} = [232.4  175.3  127.6   87.5  53.5  24.7  0]';                   %0 kPa force, N
x10{2} = [-.03   -.02   -.01   0    .01   .02]';
z10{2} = [321.3 213.3  137.1  81.7  40.5  9.3]';                                                                 %100 kPa force, N
x10{3} = [-.03   -.02   -.01     0    .01    .02   .03   .04   .05   .06]';
z10{3} = [409.8  298.2  219.6  162.2  119.1  85.8  59.7  38.8  21.9  7.9]';                                      %200 kPa force, N
x10{4} = [-.03   -.02   -.01       0    .01  .02   .03   .04   .05   .08   .11]';
z10{4} = [497.8  382.3  301.1  241.7  196.6  161.3 133  109.8  90.3  46.7  15.8]';                               %300 kPa force, N
x10{5} = [-.03   -.02   -.01    0  .01   .02     .03   .04    .05    .08   .11  .14  .15  .17]';
z10{5} = [585.4 465.5  381.6  320  273  235.8  205.5   180  158.2  106.7  67.1  34  23.8  4.6]';                 %400 kPa force, N
x10{6} = [-.02   -.01     0    .01    .02   .03   .04    .05    .08   .11  .14  .15  .17   .2]';
z10{6} = [547.8 460.9 397.2  348.3  309.3  277  249.5  225.5  166.6  118.9  77  64  38.7  2.7]';                 %500 kPa force, N
x10{7} = [-.02   -.01       0     .01    .02   .03    .04    .05    .08     .11    .14   .15   .17   .2   .22]';
z10{7} = [629.3  539.3  473.3   422.5  381.7 347.6  318.2  292.1  226.4  171.1  121.2  105.3  74.4  29.9  1.1]';  %600 kPa force, N
x10{8} = [-0.02   -.01      0    .01    .02   .03    .04    .05   .08   .11    .14    .15  .17  .2   .22]';
z10{8} = [645.4 554.8  488.4  437.3  396.1 361.6  331.8  305.4  238.4 181.6 130.1  113.7  81.8 35.6  5.7]';                 %620 kPa force, N
x10{9} = [-.01       0    .01   .02    .03   .04    .05   .08     .11    .14    .15    .16  .18  .2     .22 .23]';
z10{9} = [616.5  548.3  495.7  453.1 417.3 386.1  358.1  286.2  223.8  166.3  147.8  129.6   94  59.1  24.9   8]';                 %700 kPa force, N
x10{10} = [0      .01      .02   .03    .04    .05   .06    .08   .11   .15    .2   .24]';
z10{10} = [622.1  567.6  523.5   486  453.2  423.5 396.1  345.8  277  191.5  90.4  12.4]';                 %800 kPa force, N

y10 = cell(length(x10),1);
    for i = 1:length(x10)
        y10{i} = [Y2(i).*ones(length(x10{i}),1)]; 
    end


%Remove values with a lot of stretching or those above the Festo specified
%maximum force value
Festemp = cell2mat([x10, y10, z10]);
Fest10 = Festemp((Festemp(:,1)>=Emin10 &Festemp(:,3)<=z10maxFesto),:);

x_10 = Fest10(:,1); y_10=Fest10(:,2); z_10=Fest10(:,3);     %convert to column vectors again

%normalize
x10norm = x_10/E10max;
y10norm = y_10/P10max;
z10norm = z_10/z10max;

[XX10, YY10, ZZ10] = prepareSurfaceData(x10norm, y10norm, z10norm);

%% From Curve Fit app generate code
% For Hunt's 10mm Data:
RelativeStrain = linspace(0,1,30);
P = linspace(0,620,19);
load ForceStrainForFit.mat z



%% Load fit results from bpaFits.m (Testing_Data folder)
load bpaFitsResult.mat fitresult gof output valid XX YY ZZ

%% Create lookup tables
f_10 = fitresult{3};
[X2h,Y2h] = meshgrid(X2./E10max,Y2./Pmax);
FestoLookup10 = f_10(X2h,Y2h);

f20 = fitresult{1};
[X2g, Y1g] = meshgrid(X2./E20max,Y1./Pmax);
FestoLookup20 = f20(X2g,Y1g);

f40 = fitresult{4};
[X1h,Y1h] = meshgrid(X1./E40max,Y1./Pmax);
FestoLookup40 = f40(X1h,Y1h);

Xgr = {X2h, X2g, X1h};
Ygr = {Y2h, Y1g, Y1h};
FestoLookup = {FestoLookup10, FestoLookup20, FestoLookup40};


%% Get experimental data ready for plotting
load allData.mat Xf Yf Zf
Datmat = [Xf, Yf, Zf];
A = sortrows(Datmat,[1 2 3],'ascend');
D = cell(length(Y2),1);
buff = 1;                 %buffer around P value (in kPa) to incorporate in plot
for i = 1:length(Y2)
    D{i} = A(( round(A(:,2)*620)<=(Y2(i)+buff) & round(A(:,2)*620)>=(Y2(i)-buff)   ),:);
end

%% Create figure
Dia = ["10","20","40"];
dstr = ["model","festo","experiment"];
for k = 1: length(Dia)
    figure
    subplot(3,1,k)
    xlabel('Contraction \\epsilon^*','interpreter','tex'),ylabel('Force, F^*','interpreter','tex')
    strTit = sprintf('\\phi%s mm BPA Force-Pressure-Contraction relationship',Dia(k));
    title(strTit,'interpreter','tex')
    hold on
    for i = 1:length(x40norm)
        for j = 1: length(dstr)
            str(i,j) = sprintf('P^*=$-.3f, %s',Ygr{k}(i,j),dstr{j});
        end
        scatter(x40norm{i},z40norm{i},'o','DisplayName',str(1))
        plot(Xgr{k},FestoLookup{k}(i,:),'DisplayName',str(2))
        if k == 1
            for r = 1:length(D)
                if ~isempty(D{r})
                    x = D{r}(:,1);
                    z = D{r}(:,3);
                    scatter(x,z,[],'DisplayName',str(3))
                else
                end
            end
        else
        end
    end
    lgd(k) = legend;
    title(lgd(k),'P^* value, data source')
    hold off
end

% 
% figure
% xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Force, $N$','interpreter','latex')
% title('\bf 20 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')
% hold on
% plot(X2,FestoLookup20(8,:)*z20max,'DisplayName','620 kPa, model')
% plot(x20{8},z20{8},'o','DisplayName','620 kPa,  Festo')
% plot(X2,FestoLookup20(7,:)*z20max,'DisplayName','600 kPa, model')
% plot(x20{7},z20{7},'o','DisplayName','600 kPa,  Festo')
% plot(X2,FestoLookup20(6,:)*z20max,'DisplayName','500 kPa, model')
% plot(x20{6},z20{6},[],'o','DisplayName','500 kPa,  Festo')
% plot(X2,FestoLookup20(5,:)*z20max,'DisplayName','400 kPa, model')
% plot(x20{5},z20{5},'o','DisplayName','400 kPa,  Festo')
% plot(X2,FestoLookup20(4,:)*z20max,'DisplayName','300 kPa, model')
% plot(x20{4},z20{4},'o','DisplayName','300 kPa,  Festo')
% plot(X2,FestoLookup20(3,:)*z20max,'DisplayName','200 kPa, model')
% plot(x20{3},z20{3},'o','DisplayName','200 kPa,  Festo')
% plot(X2,FestoLookup20(2,:)*z20max,'DisplayName','100 kPa, model')
% plot(x20{2},z20{2},'o','DisplayName','100 kPa,  Festo')
% plot(X2,FestoLookup20(1,:)*z20max,'DisplayName','    0 kPa, model')
% plot(x20{1},z20{1},'o','DisplayName','    0 kPa,  Festo')
% lgd20 = legend;
% title(lgd20,'\bf Pressure')
% hold off

% figure
% surf(X2,Y1,FestoLookup20*z20max)
% xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Pressure, kPA','interpreter','latex'),zlabel('\bf Force, N','interpreter','latex')
% title('\bf 20 $ mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')
% 
% figure
% surf(X1,Y1,FestoLookup40)
% xlabel('\bf Contraction','interpreter','latex'),ylabel('\bf Pressure, $kPA$','interpreter','latex'),zlabel('\bf Force, $N$','interpreter','latex')
% title('40 $mm$ BPA Force-Pressure-Contraction relationship','interpreter','latex')
%
