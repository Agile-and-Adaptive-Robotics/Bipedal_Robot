%% Import data from Lawrence's two sets of tests, compare with Ben's. Upload for curve fit
clear;
clc;
close all

load FestoData.mat
%% Lawrence's first set of tests
load allData_Ben.mat Cut li lo l620 maxE relE T test Q
Cuts1 = Cut;      %cut lengths
Ts1 = T;          %data
lis1 = li;        %contracted length
lmin1 = l620;     %minimum length
restL1 = lo;      %resting lengths
maxE1 = maxE;     %maximum strain
relE1 = relE;     %relative strain
tests1 = test;    %cell of strings indicating test number
Qs1 = Q;          %full set of data from all the tests

clear Cut li lo l620 maxE relE T test Q

%% Lawrence's first set of tests
load allData_Nov22_Ben.mat Cut li lo l620 maxE relE T vars Q
Cuts2 = Cut;      %cut lengths
Ts2 = T;          %data
lis2 = li;        %contracted length
lmin2 = l620;     %minimum length
restL2 = lo;      %resting lengths
maxE2 = maxE;     %maximum strain
relE2 = relE;     %relative strain
tests2 = vars;    %cell of strings indicating test number
Qs2 = Q;          %full set of data from all the tests

clear Cut li lo l620 maxE relE T vars Q

%% Combine both into cells
Cut = {Cuts1, Cuts2};
T = {Ts1, Ts2};
li = {lis1, lis2};
lmin = {lis1, lis2};
restL = {restL1, restL2};
maxE = {maxE1, maxE2};
relE = {relE1, relE2};
tests = {tests1, tests2};
Q = {Qs1, Qs2};


%% Find max. forces
Ymax = [100, 200, 300, 400, 500, 620];
meth = 'linear';
extrap = 'extrap';

maxP = 800;
minF = 5;


maxF = cell(1,length(li));
for a = 1:length(li)
   maxF{a} = cell(length(Cut{a}), length(tests{a}));
    for i = 1:length(Cut{a})
        for k = 1:length(tests{a})
        M = T{a}{i,1,k};       %first column, looking at resting length force
        M = M( (M(:,1)==0 &M(:,3)>minF &M(:,2)<maxP),1:3);
        pres = M(:,2);
        fz = M(:,3);
        maxF{a}{i,k} = interp1(pres,fz,Ymax,meth,extrap);
        end
       gg = vertcat(maxF{a}{i,:});
       Fmax{i,a} = mean(gg,1,'omitnan');
    end

end

B = vertcat(Fmax{:});
F = max(B,2);    %I think this doesn't do anything?
F_max = max(F,[],2);

% [Fmax13, Fmax23, Fmax27, Fmax29,Fmax30,Fmax10,Fmax15,Fmax20,Fmax25,Fmax30_2,Fmax40,Fmax45_2,Fmax52_2]'=F(:)';
keep = [3 1]; %keep test 9 (keep(2)) from first set of data, test 8 from second set (keep(1))
Ts = vertcat(T{1}(:,:,keep(1)),T{2}(:,:,keep(2)));  
[a,b] = size(Ts);
V = cell(size(Ts));
Pmax = 620;
for i = 2:a             %remove data for first cut length (l_rest = 120 mm)
    for j = 1:b
            if ~isempty(Ts{i,j})
            V{i,j} = Ts{i,j}((Ts{i,j}(:,3)>minF &Ts{i,j}(:,2)<maxP),1:3);    
            V{i,j} = [Ts{i,j}(:,1), Ts{i,j}(:,2)./Pmax, Ts{i,j}(:,3)./F_max(i)];
            else
            end
    end
end

%Repeat above but for full size data
Qs = vertcat(Q{1}(:,:,keep(1)),Q{2}(:,:,keep(2)));  
U = cell(size(Qs));
for i = 2:a            %remove data for first cut length (l_rest = 120 mm)
    for j = 1:b
            if ~isempty(Qs{i,j})
            U{i,j} = Qs{i,j}((Qs{i,j}(:,3)>minF &Qs{i,j}(:,2)<maxP),1:3);    
            U{i,j} = [Qs{i,j}(:,1), Qs{i,j}(:,2)./Pmax, Qs{i,j}(:,3)./F_max(i)];
            else
            end
    end
end

H = [1 1 0;...    %full strain, full pressure, no force
    0 1 1;...     %no strain, full pressure, full force
    0 0 0];       %no strain, no pressure, no force

%Reduced data
R = vertcat(V{:});
R = [R; H];
R1 = vertcat(V{1:length(Cuts1),:});
R1 = [R1; H];
R2 = vertcat(V{(length(Cuts1)+1):end,:});
R2 = [R2; H];
Rx = R(:,1); Ry = R(:,2); Rz = R(:,3);
R1x = R1(:,1); R1y = R1(:,2); R1z = R1(:,3);
R2x = R2(:,1); R2y = R2(:,2); R2z = R2(:,3);

%Full size data
S = vertcat(U{:});
S = [S; H];
S1 = vertcat(U{1:length(Cuts1),:});
S1 = [S1; H];
S2 = vertcat(U{(length(Cuts1)+1):end,:});
S2 = [S2; H];
Sx = S(:,1); Sy = S(:,2); Sz = S(:,3);
S1x = S1(:,1); S1y = S1(:,2); S1z = S1(:,3);
S2x = S2(:,1); S2y = S2(:,2); S2z = S2(:,3);
%% Add Ben's data
%Fmax from experiments, using gridded or scattered interpolation
Fmax112 = 350.86;
Fmax415 = 444.82;
Fmax455 = 440.74;
Fmax490 = 451.15;
Fmax518 = 456.17;

rawdata11cm = [325.164999	620	0.008928571	0.055555556;
                252.6589869	620	0.044642857	0.277777778;
                135.6707588	620	0.089285714	0.555555556;
                89.40925416	620	0.125	0.777777778;
                13.3446648	620	0.169642857	1.055555556
                        0	620	0.160714286	1];
% gg = sortrows([rawdata11cm(:,4), rawdata11cm(:,1)]);
% max112 = griddedInterpolant(gg(:,1),gg(:,2),'linear','linear');
% Fmax112 = max112(0);
% disp(Fmax112)
data11cm = [rawdata11cm(:,4), rawdata11cm(:,2), rawdata11cm(:,1)/Fmax112
            0, 0, 0];

% zcomp = max112([0; gg(:,1)]);
% figure
% hold on
% scatter(rawdata11cm(:,4),rawdata11cm(:,1),[],'bo');
% plot([0; gg(:,1)],zcomp, '-.r');
% hold off

rawdata42cm = [444.82216	620	0.002409639	0.014492754;
                358.9714831	620	0.019277108	0.115942029;
                275.7897392	620	0.031325301	0.188405797;
                200.169972	620	0.06746988	0.405797101;
                103.1987411	620	0.113253012	0.68115942;
                40.92363872	620	0.151807229	0.913043478;
                       0	620	0.16626506	1];
% gg = sortrows([rawdata42cm(:,4),rawdata42cm(:,1)]);
% max415 = griddedInterpolant(gg(:,1),gg(:,2),'linear','makima');
% Fmax415 = max415(0);
data42cm = [rawdata42cm(:,4), rawdata42cm(:,2), rawdata42cm(:,1)/Fmax415
            0, 0, 0];

% zcomp = max415([0; gg(:,1)]);
% figure
% hold on
% scatter(rawdata42cm(:,4),rawdata42cm(:,1),[],'bo');
% plot([0; gg(:,1)],zcomp, '-.r');
% hold off

rawdata45cm = [26	620	0.161147903	1
61	620	0.14790287	0.917808219
150	620	0.097130243	0.602739726
96	620	0.135761589	0.842465753
100	620	0.130242826	0.808219178
119	620	0.125827815	0.780821918
158	620	0.105960265	0.657534247
217.0732141	620	0.079470199	0.493150685
262.0002522	620	0.057395143	0.356164384
335.3959086	620	0.033112583	0.205479452
407.0122764	620	0.017660044	0.109589041
435.9257168	620	0.002207506	0.01369863
404.7881656	620	0.006622517	0.04109589
408.3467429	620	0.004415011	0.02739726
335.3959086	620	0.017660044	0.109589041
257.9968528	620	0.041942605	0.260273973
203.2837271	620	0.070640177	0.438356164
0	620	0.161147903	1];
% gg = sortrows([rawdata45cm(:,4),rawdata45cm(:,1)],1);
% [p, gof, output] = fit(gg(:,1),gg(:,2),'poly3','Normalize','on','Robust','on');
% ggRsz = imresize(gg, [length(gg)/2, 2]);
% gg = groupsummary(gg, gg(:,1), 'mean');
% max455 = griddedInterpolant(ggRsz(:,1),ggRsz(:,2),'linear','linear');
% Fmax455 = p(0);
data45cm = [rawdata45cm(:,4), rawdata45cm(:,2), rawdata45cm(:,1)/Fmax455
            0, 0, 0];

% zcomp = max455([0; gg(:,1)]);
% xi = 0:0.1:1;
% zi = feval(p,xi);
% figure
% hold on
% scatter(rawdata45cm(:,4),rawdata45cm(:,1),5,'bo');
% plot([0; gg(:,1)],zcomp, '-.r');
% scatter(ggRsz(:,1),ggRsz(:,2),10,'d','b');
% plot(xi,zi, '--k');
% hold off


rawdata49cm = [464	620	-0.002040816	-0.010869565
416	616	0.008163265	0.043478261
308.77	618	0.030612245	0.163043478
356	620	0.02244898	0.119565217
300	618	0.03877551	0.206521739
268	617	0.048979592	0.260869565
241	617	0.06122449	0.326086957
215	617	0.069387755	0.369565217
194	617	0.081632653	0.434782609
165	617	0.095918367	0.510869565
144	617	0.106122449	0.565217391
107	617	0.128571429	0.684782609
97	617	0.132653061	0.706521739
75	617	0.146938776	0.782608696
58	617	0.157142857	0.836956522
45	617	0.165306122	0.880434783
20	617	0.17755102	0.945652174
2.29	617	0.185714286	0.989130435
0   620	0.187755102	1
0 0 0 0];
gg = sortrows([rawdata49cm(:,4),rawdata49cm(:,2),rawdata49cm(:,1)],1);
max490 = scatteredInterpolant(gg(:,1), gg(:,2), gg(:,3),'linear','linear');
Fmax490 = max490(0,620);
data49cm = [rawdata49cm(:,4), rawdata49cm(:,2), rawdata49cm(:,1)/Fmax490];

xi = [0; gg(:,1)];
yi = 620*ones(length(xi),1);
zcomp = max490(xi, yi);
figure
hold on
scatter3(gg(:,1),gg(:,2),gg(:,3),[],'bo');
plot3(xi,yi,zcomp, '-.r');
hold off

rawdata52cm = [0	619	0.166023166	1
% -6	618	0.167953668	1.011627907
% -5.3	618	0.171814672	1.034883721
60	616	0.150579151	0.906976744
427.6	615	0.005791506	0.034883721
347	616.8	0.023166023	0.139534884
296	618	0.038610039	0.23255814
245	617	0.055984556	0.337209302
191	617	0.079150579	0.476744186
55	617	0.152509653	0.918604651
99	617	0.135135135	0.813953488
132	617	0.115830116	0.697674419
160	617	0.102316602	0.61627907
195	617	0.086872587	0.523255814
450	617	0.003861004	0.023255814
420	617	0.011583012	0.069767442
328	617	0.030888031	0.186046512
378	617	0.019305019	0.11627907
230	617	0.063706564	0.38372093
278	617	0.05019305	0.302325581
429	617	0	0
0	0	0	0];

gg = sortrows([rawdata52cm(:,4), rawdata52cm(:,2), rawdata52cm(:,1)],1);
[gx, gy, gz] = prepareSurfaceData( gg(:,1), gg(:,2)/620, gg(:,3) );

% Set up fittype and options.
ft = fittype( 'a*x^3+b*x^2+c*x+f*y', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-800 1400 -1115 492];
opts.Robust = 'Bisquare';

% Fit model to data.
[max518, gof] = fit( [gx, gy], gz, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 2' );
h = plot( max518, [gx, gy], gz );
legend( h, 'untitled fit 2', 'gz vs. gx, gy', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'gx', 'Interpreter', 'none' );
ylabel( 'gy', 'Interpreter', 'none' );
zlabel( 'gz', 'Interpreter', 'none' );
grid on
view( -51.0, 38.5 );

Fmax518 = max518(0,1);
data52cm = [rawdata52cm(:,4), rawdata52cm(:,2), rawdata52cm(:,1)/Fmax518];

xi = [0; gg(:,1)];
yi = ones(length(xi),1);
zcomp = max518(xi, yi);
figure
hold on
scatter3(gg(:,1),gg(:,2)/620,gg(:,3),[],'bo');
plot3(xi,yi,zcomp, '-.r');
hold off

%% Combine all data and do a 3d scatter plot. 
%addnoise to anchor surface at three points
sz = 15;                   % size of anchor surface
err = 0.001;                 % ammount of error
r110 = [1-err + 2*err * rand(sz,2), -err + 2*err * rand(sz,1)].*[1 1 1];
r011 = [-err + 2*err * rand(sz,1), 1-err + 2*err * rand(sz,2)].*[1 1 1];
r000 = -0.001 + 0.002 * rand(sz,3);
anchor = [r110; r011; r000];
ancX = anchor(:,1); ancY = anchor(:,2); ancZ = anchor(:,3);


% Ben_data = [data11cm; data42cm; data45cm; data49cm; data52cm];
Ben_data = [data42cm; data49cm; data52cm];
Ben_data = [Ben_data(:,1), Ben_data(:,2)./620, Ben_data(:,3)];
BenX = Ben_data(:,1); BenY = Ben_data(:,2); BenZ = Ben_data(:,3);

allData = [ R1;...
%             R2;...
            Ben_data;...
            anchor];


w = [0.5*ones(length(R1),1); 
%     0.66*ones(length(R2),1); 
    2*ones(length(Ben_data),1); 
    ones(length(anchor),1)]; %set weights

[XX, YY, ZZ, WW] = prepareSurfaceData(allData(:,1), allData(:,2),allData(:,3),w);
% XX = allData(:,1); YY=allData(:,2); ZZ=allData(:,3);
%Ypres = YY*620;   %YY repressurized


w = [0.5*ones(length(S1),1); 
%     0.66*ones(length(S2),1); 
    2*ones(length(Ben_data),1); 
    ones(length(anchor),1)]; %set weights

bigData = [ S1;...
%             S2;...
            Ben_data;...
            anchor];
        
[Xf, Yf, Zf, Wf] = prepareSurfaceData(bigData(:,1), bigData(:,2), bigData(:,3),w);

figure
hold on
scatter3(Xf, Yf, Zf, 'DisplayName', 'All data, full sized')
scatter3(XX, YY, ZZ,'filled', 'DisplayName', 'All data, resized')
hold off
xlabel('Relative Strain')
ylabel('Pressure (normalized)')
zlabel('Force (normalized)')
title('10mm force, normalized')
lgd1 = legend;

figure
hold on
scatter3(R1(:,1), R1(:,2), R1(:,3),'DisplayName','Lawrence 1')
scatter3(R2(:,1), R2(:,2), R2(:,3),'DisplayName','Lawrence 2')
scatter3(Ben_data(:,1),Ben_data(:,2),Ben_data(:,3),'DisplayName','Ben data')
scatter3(anchor(:,1), anchor(:,2), anchor(:,3),'DisplayName','anchors')
hold off
xlabel('Relative Strain')
ylabel('Pressure (normalized)')
zlabel('Force (normalized)')
title('10mm force, normalized')
lgd2 = legend;

% Assuming x, y, z are column vectors of the same length
n = length(XX);
rng('default'); % For reproducibility
idx = randperm(n);

% 80% training, 20% validation
numTrain = round(0.8 * n);
trainIdx = idx(1:numTrain);
valIdx = idx(numTrain+1:end);

% Training data
xTrain = XX(trainIdx);
yTrain = YY(trainIdx);
zTrain = ZZ(trainIdx);

% Validation data
xVal = XX(valIdx);
yVal = YY(valIdx);
zVal = ZZ(valIdx);

save allData.mat    XX YY ZZ WW... 
                    Ben_data BenX BenY BenZ...
                    R Rx Ry Rz...
                    R1 R1x R1y R1z...
                    R2 R2x R2y R2z...
                    anchor ancX ancY ancZ...
                    Xf Yf Zf Wf...
                    F F_max...
                    S Sx Sy Sz...
                    S1 S1x S1y S1z...
                    S2 S2x S2y S2z...
                    xTrain yTrain zTrain...
                    xVal yVal zVal