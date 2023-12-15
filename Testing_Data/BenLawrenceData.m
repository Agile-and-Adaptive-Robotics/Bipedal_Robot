%% Import data from Lawrence's two sets of tests, compare with Ben's. Upload for curve fit
clear; clc; close all

%% Lawrence's first set of tests
load allData_Ben.mat Cut li lo l620 maxE relE T test
Cuts1 = Cut;       %cut lengths
Ts1 = T;           %data
lis1 = li;        %contracted length
lmin1 = l620;     %minimum length
restL1 = lo;      %resting lengths
maxE1 = maxE;     %maximum strain
relE1 = relE;     %relative strain
tests1 = test;  %cell of strings indicating test number

clear Cut li lo l620 maxE relE T test

%% Lawrence's first set of tests
load allData_Nov22_Ben.mat Cut li lo l620 maxE relE T vars
Cuts2 = Cut;       %cut lengths
Ts2 = T;           %data
lis2 = li;        %contracted length
lmin2 = l620;     %minimum length
restL2 = lo;      %resting lengths
maxE2 = maxE;     %maximum strain
relE2 = relE;     %relative strain
tests2 = vars;  %cell of strings indicating test number

clear Cut li lo l620 maxE relE T vars

%% Combine both into cells
Cut = {Cuts1, Cuts2};
T = {Ts1, Ts2};
li = {lis1, lis2};
lmin = {lis1, lis2};
restL = {restL1, restL2};
maxE = {maxE1, maxE2};
relE = {relE1, relE2};
tests = {tests1, tests2};


%% Find max. forces
Ymax = [100, 200, 300, 400, 500, 620];
meth = 'linear';
extrap = 'extrap';

maxP = 800;
minF = 20;


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
keep = [2 1]; %keep test 9 from first set of data, test 8 from second set
Ts = vertcat(T{1}(:,:,keep(1)),T{2}(:,:,keep(2)));  
[a,b] = size(Ts);
V = cell(size(Ts));
Pmax = 620;
for i = 1:a
    for j = 1:b
            if ~isempty(Ts{i,j})
            V{i,j} = Ts{i,j}((Ts{i,j}(:,3)>minF &Ts{i,j}(:,2)<maxP),1:3);    
            V{i,j} = [Ts{i,j}(:,1), Ts{i,j}(:,2)./Pmax, Ts{i,j}(:,3)./F_max(i)];
            else
            end
    end
end

R = vertcat(V{:});
R1 = vertcat(V{1:length(Cuts1),:});
R2 = vertcat(V{(length(Cuts1)+1):end,:});

Rx = R(:,1); Ry = R(:,2); Rz = R(:,3);
R1x = R1(:,1); R1y = R1(:,2); R1z = R1(:,3);
R2x = R2(:,1); R2y = R2(:,2); R2z = R2(:,3);
%% Add Ben's data
%Fmax from experiments, using interpolation
Fmax112 = 334.74;
Fmax415 = 445;
Fmax455 = 460;
Fmax490 = 458.6;
Fmax518 = 455.83;

rawdata11cm = [325.164999	620	0.008928571	0.055555556;
                252.6589869	620	0.044642857	0.277777778;
                135.6707588	620	0.089285714	0.555555556;
                89.40925416	620	0.125	0.777777778;
                13.3446648	620	0.169642857	1.055555556];
data11cm = [rawdata11cm(:,4), rawdata11cm(:,2), rawdata11cm(:,1)/Fmax112];

rawdata42cm = [444.82216	620	0.002409639	0.014492754;
                358.9714831	620	0.019277108	0.115942029;
                275.7897392	620	0.031325301	0.188405797;
                200.169972	620	0.06746988	0.405797101;
                103.1987411	620	0.113253012	0.68115942;
                    40.92363872	620	0.151807229	0.913043478;
                                0	620	0.16626506	1];
data42cm = [rawdata42cm(:,4), rawdata42cm(:,2), rawdata42cm(:,1)/Fmax415];

rawdata45cm = [26	620	0.164835165	1.041666667
61	620	0.151648352	0.958333333
150	620	0.101098901	0.638888889
96	620	0.13956044	0.881944444
100	620	0.134065934	0.847222222
119	620	0.12967033	0.819444444
158	620	0.10989011	0.694444444
217.0732141	620	0.083516484	0.527777778
262.0002522	620	0.061538462	0.388888889
335.3959086	620	0.037362637	0.236111111
407.0122764	620	0.021978022	0.138888889
435.9257168	620	0.006593407	0.041666667
404.7881656	620	0.010989011	0.069444444
408.3467429	620	0.008791209	0.055555556
335.3959086	620	0.021978022	0.138888889
257.9968528	620	0.046153846	0.291666667
203.2837271	620	0.074725275	0.472222222];
data45cm = [rawdata45cm(:,4), rawdata45cm(:,2), rawdata45cm(:,1)/Fmax455];

rawdata49cm = [12	618	0.173469388	0.923913043
1.11	618	0.175510204	0.934782609
464	620	-0.002040816	-0.010869565
433	616	0.008163265	0.043478261
320	618	0.030612245	0.163043478
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
0   620	0.187755102	1];
data49cm = [rawdata49cm(:,4), rawdata49cm(:,2), rawdata49cm(:,1)/Fmax490];

rawdata52cm = [0	619	0.166023166	1
6	618	0.167953668	1.011627907
-5.3	618	0.171814672	1.034883721
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
429	617	0	0];
data52cm = [rawdata52cm(:,4), rawdata52cm(:,2), rawdata52cm(:,1)/Fmax518];

%% Combine all data and do a 3d scatter plot. 
%addnoise to anchor surface at three points
sz = 15;                   % size of anchor surface
err = 0.001;                 % ammount of error
r110 = [1-err + 2*err * rand(sz,2), -err + 2*err * rand(sz,1)].*[1 1 1];
r011 = [-err + 2*err * rand(sz,1), 1-err + 2*err * rand(sz,2)].*[1 1 1];
r000 = -0.01 + 0.02 * rand(sz,3);
anchor = [r110; r011; r000];
ancX = anchor(:,1); ancY = anchor(:,2); ancZ = anchor(:,3);


Ben_data = [data11cm; data42cm; data45cm; data49cm; data52cm];
Ben_data = [Ben_data(:,1), Ben_data(:,2)./620, Ben_data(:,3)];
BenX = Ben_data(:,1); BenY = Ben_data(:,2); BenZ = Ben_data(:,3);

allData = [ R1;...
            R2;...
            Ben_data;...
            anchor];


w = [0.33*ones(length(R1),1); 
    0.66*ones(length(R2),1); 
    2*ones(length(Ben_data),1); 
    0.25*ones(length(anchor),1)]; %set weights

[XX, YY, ZZ, WW] = prepareSurfaceData(allData(:,1), allData(:,2),allData(:,3),w);
% XX = allData(:,1); YY=allData(:,2); ZZ=allData(:,3);
%Ypres = YY*620;   %YY repressurized

figure
hold on
scatter3(XX, YY, ZZ)
hold off
xlabel('Relative Strain')
ylabel('Pressure (normalized)')
zlabel('Force (normalized)')
title('10mm force, normalized')

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
lgd = legend;

save allData.mat    XX YY ZZ... 
                    Ben_data BenX BenY BenZ...
                    R Rx Ry Rz...
                    R1 R1x R1y R1z...
                    R2 R2x R2y R2z...
                    anchor ancX ancY ancZ...
                    w...
                    F F_max