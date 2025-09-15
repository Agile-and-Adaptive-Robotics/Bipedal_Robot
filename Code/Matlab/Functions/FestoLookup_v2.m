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
% E10min= -0.03;
% E20min= -0.04;
% E40min= -0.05;

%disregard sretching more than 2% resting length
E10min= -0.02;
E20min= -0.02;
E40min= -0.02;

X1 = linspace(-0.02,0.32,35);   %Strain range for interpolation, 40mm
X2 = linspace(-0.02,0.28,31);   %Strain range for interpolation, 10 & 20mm
%% 40mm dia BPA
x40 = cell(length(Y1),1);
z40 = cell(length(Y1),1);

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


y40 = cell(length(Y1),1);
    for i = 1:length(x40)
        y40{i} = [Y1(i).*ones(length(x40{i}),1)]; 
    end

%Remove values with a lot of stretching or those above the Festo specified
%value
Festemp = cell2mat([x40, y40, z40]);
Fest40 = Festemp((Festemp(:,1)>=E40min &Festemp(:,3)<=z40maxFesto &Festemp(:,3)>=0),:);

x_40 = Fest40(:,1); y_40=Fest40(:,2); z_40=Fest40(:,3);     %convert to column vectors again


x40norm = x_40/E40max;
y40norm = y_40/Pmax;
z40norm = z_40/z40max;

[XX40, YY40, ZZ40] = prepareSurfaceData(x40norm, y40norm, z40norm);
%% 20 mm BPA
x20 = cell(length(Y1),1);
z20 = cell(length(Y1),1);

x20{1} = [-.04 -.03 -.02 -.01 0 0.03 .125 .25]';
z20{1} = [ 684  416  229   96 0 -160 -270 -287]';                   %0 kPa force, N
x20{2} = [-.04 -.03 -.02 -.01 0   0.02 0.06 0.17 .25]';
z20{2} = [ 940  672  486  352 253 123    -8 -147 -218]';            %100 kPa force, N
x20{3} = [-.04 -.02 -.01 0.01 .05 .12 .155 .25]';
z20{3} = [1194  743  608  429 240 112 8.7  -150]';                  %200 kPa force, N
x20{4} = [-.04 -.03 -.01 .01 .04 .12 .215 .25]';
z20{4} = [1451 1187  864 678 507 241  0.5 -84]';                    %300 kPa force, N
x20{5} = [-.04 -.02    0 .01 .04 .09 .16 .24 .25]';
z20{5} = [1707 1257 1012 926 737 518 271  13 -19]';                 %400 kPa force, N
x20{6} = [-.04 -.019    0 0.01 .04 .08 .13 .195 .22 .25  .26]';
z20{6} = [1963  1498 1264 1173 966 756 531  262 162  45  5.8]';     %500 kPa force, N
x20{7} = [-.04 0.00 .002 .035 .065 .11 .14 .17  .20  .22 .25 .273]';
z20{7} = [2219 1515 1495 1228 1040 791 637 487  342  247 107  0.9]';       %600 kPa force, N
x20{8} = [-.04 0.00 .007 .035 .065 .11  .14 .17  .20  .22  .25 .275]';
z20{8} = [2270 1565 1497 1274 1082 826  666 512  362  264  119  0.3]';     %620 kPa force, N


y20 = cell(length(x20),1);
    for i = 1:length(x20)
        y20{i} = [Y1(i).*ones(length(x20{i}),1)]; 
    end

%Remove values with a lot of stretching or those above the Festo specified
%value
Festemp = cell2mat([x20, y20, z20]);
Fest20 = Festemp((Festemp(:,1)>=E20min &Festemp(:,3)<=z20maxFesto &Festemp(:,3)>=0),:);

x_20 = Fest20(:,1); y_20=Fest20(:,2); z_20=Fest20(:,3);     %convert to column vectors again

x20norm = x_20/E20max;
y20norm = y_20/Pmax;
z20norm = z_20/z20max;

[XX20, YY20, ZZ20] = prepareSurfaceData(x20norm, y20norm, z20norm);
%% 10 mm BPA
X3 = linspace(-0.03,0.25,29);   %Strain range for interpolation
x10 = cell(length(Y2),1);
z10 = cell(length(Y2),1);

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

y10 = cell(length(Y2),1);
    for i = 1:length(x10)
        y10{i} = [Y2(i).*ones(length(x10{i}),1)]; 
    end


%Remove values with a lot of stretching or those above the Festo specified
%maximum force value
Festemp = cell2mat([x10, y10, z10]);
Fest10 = Festemp((Festemp(:,1)>=E10min &Festemp(:,3)<=z10maxFesto),:);

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
load bpaFitsResult.mat fitresult gof output valid XX YY ZZ Ax Ay Az

%% Create lookup tables
Yp = [200, 300, 400, 500, 620]; %"Y prime", Remove 600 kPA

f_10 = fitresult{1};                            %10 mm, using data
[X3g,Y2g] = meshgrid(X3./E10max,Yp./P10max);
FestoLookup10 = f_10(X3g,Y2g);

f20 = fitresult{2};                             %20 mm, using data
[X2g, Y1g] = meshgrid(X2./E20max,Yp./Pmax);
FestoLookup20 = f20(X2g,Y1g);

f40 = fitresult{3};                             %40 mm, using data
[X1g,Y1h] = meshgrid(X1./E40max,Yp./Pmax);
FestoLookup40 = f40(X1g,Y1h);

Xgr = {X3g, X2g, X1g};
Ygr = {Y2g, Y1g, Y1h};
FestoLookup = {FestoLookup10, FestoLookup20, FestoLookup40};

% save FestoLookup.mat f_10 f20 f40


%% Get 10 mm experimental data ready for plotting
load allData.mat Xf Yf Zf
Datmat = [Xf, Yf, Zf];
A = sortrows(Datmat,[1 2 3],'ascend');
D = cell(length(Yp),1);
buff = 1;                 %buffer around P value (in kPa) to incorporate in plot
for i = 1:length(Yp)
    D{i} = A(( round(A(:,2)*620)<=(Yp(i)+buff) & round(A(:,2)*620)>=(Yp(i)-buff)   ),:);
end

%% Get 20 mm experimental data ready for plotting
load data20mm_sorted.mat Ax Ay Az
Datmat = [Ax, Ay, Az];
A = sortrows(Datmat,[1 2 3],'ascend');
E = cell(length(Yp),1);
buff = 5;                 %buffer around P value (in kPa) to incorporate in plot
for i = 1:length(Yp)
    E{i} = A(( round(A(:,2)*620)<=(Yp(i)+buff) & round(A(:,2)*620)>=(Yp(i)-buff)   ),:);
end

%% "Get rid of 0, 100, 600, 700, and 800 kPa Festo data, keep 620" - Dr. Hunt
x10p = [x10(3:6); x10(8)];
y10p = [y10(3:6); y10(8)];
z10p = [z10(3:6); z10(8)];

x20p = [x20(3:6); x20(8)];
y20p = [y20(3:6); y20(8)];
z20p = [z20(3:6); z20(8)];

x40p = [x40(3:6); x40(8)];
y40p = [y40(3:6); y40(8)];
z40p = [z40(3:6); z40(8)];

%% Addition figure script setup
Dia = ["10","20","40"];
dstr = ["festo","model","experiment"];

Xc = {x10p, x20p, x40p};
Yc = {y10p, y20p, y40p};
Zc = {z10p, z20p, z40p};

Emax = [E10max, E20max, E40max];
Emin = [E10min, E20min, E40min];
Pmax = [P10max, Pmax, Pmax];
Zmax = [z10max, z20max, z40max];


%% Figure (for 10 mm & 20 mm)

    figure
    tileLabels = {'(A)', '(B)', '(C)', '(D)'};
    % Annotation positions [x, y] in normalized figure units
    xAnn = [0, 0];
    yAnn = [0.94, 0.45];
    for k = 1: 2
    a(k) = subplot(2,1,k);
    xlabel('Contraction, \epsilon^*','interpreter','tex'),ylabel('Force, F^*','interpreter','tex')
    strTit = sprintf('\\phi%s mm BPA Force-Pressure-Contraction relationship',Dia(k));
    title(strTit,'interpreter','tex')
    hold on
    for i = 1:length(Xc{k})
        for j = 1: length(dstr)
            str{i,j} = sprintf('P=%.0f kPa, %s',Ygr{k}(i)*620,dstr{j});
        end
        sc{i,k} = plot(Xc{k}{i}/Emax(k),Zc{k}{i}/Zmax(k),'--', 'DisplayName',str{i,1});
        sc{i,k}.SeriesIndex = i;
        pl{i,k} = plot(Xgr{k}(i,:),FestoLookup{k}(i,:),'DisplayName',str{i,2});
        pl{i,k}.SeriesIndex = i;
        if k == 1
            if ~isempty(D{i})
                x = D{i}(:,1);
                z = D{i}(:,3);
                sk{i,k} = scatter(x,z,120,'.','DisplayName',str{i,3});
            else
                sk{i,k} = scatter([],[],'.','DisplayName','');
            end
            sk{i,k}.SeriesIndex = i; 
        elseif k==2
            if ~isempty(E{i})
                x = E{i}(:,1);
                z = E{i}(:,3);
                sk{i,k} = scatter(x,z,120,'.','DisplayName',str{i,3});
            else
                sk{i,k} = scatter([],[],'.','DisplayName','');
            end
            sk{i,k}.SeriesIndex = i; 
        else
        end
    end
    a(k).XLim = [0 1];
    a(k).YLim = [0 max(max(Ygr{2}))];
    set(gca, 'FontSize', 10, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    lgd(k) = legend;
    lgd(k).Location = 'northeast';
    lgd(k).Orientation = 'horizontal';
    lgd(k).NumColumns = 3;
    title(lgd(k),'P value and data source')

    hold off
    
    end
    
    for j = 1:2
    annotation(gcf, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 10, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    end
    


