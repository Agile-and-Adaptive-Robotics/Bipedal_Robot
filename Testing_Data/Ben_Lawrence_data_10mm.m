%% After running LoadAll10mmData from Lawrence's muscle sensory branch
clear data13cm data23cm data27cm data 29cm data30cm

%% Plot tests 9 and 10, only Pressurizing
 testnumber = 9; %choosing test numner 9 as the standard
 testnum = 10;   %adding test 10 b/c Lawrence says it's good
%Mean maximum force at 620kPa, from tests 9 and 10, respecitively.  
Fmax13 = 341.2402; %mean([341.9415  340.5389]);
Fmax23 = 377.8612; %mean([378.1035 377.6189]);
Fmax27 = 383.3393; %mean([383.1724 383.5062]);
Fmax29 = 419.1045; %mean([419.1600 419.0490]);
Fmax30 = 406.9386; %mean([407.2156 406.6616]);
maxP = 800;
minF = 0;

%13cm
figure %11
subplot 321
hold on
for a = 1:length(vals_13cm)
    b=vals_13cm(a);
    c = [0 0.2 0.4 0.6];
    data13cm_t9 = AllBPA10mm13cm(AllBPA10mm13cm(:,5)==13& AllBPA10mm13cm(:,6)==b&AllBPA10mm13cm(:,7)==testnumber&AllBPA10mm13cm(:,8)==1,:);
    data13cm_t10 = AllBPA10mm13cm(AllBPA10mm13cm(:,5)==13& AllBPA10mm13cm(:,6)==b&AllBPA10mm13cm(:,7)==testnum&AllBPA10mm13cm(:,8)==1,:);
    txt = sprintf('%d',c(a));
    plot(data13cm_t9(:,2),data13cm_t9(:,1)/Fmax13,'DisplayName',txt)
    plot(data13cm_t10(:,2),data13cm_t10(:,1)/Fmax13,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 13cm all Kinks(Test 9&10)')
data13cm_t9 = [AllBPA10mm13cm_P(AllBPA10mm13cm_P(:,7)==testnumber,13), AllBPA10mm13cm_P(AllBPA10mm13cm_P(:,7)==testnumber,2), AllBPA10mm13cm_P(AllBPA10mm13cm_P(:,7)==testnumber,1)];                     
data13cm_t10 = [AllBPA10mm13cm_P(AllBPA10mm13cm_P(:,7)==testnum,13), AllBPA10mm13cm_P(AllBPA10mm13cm_P(:,7)==testnum,2), AllBPA10mm13cm_P(AllBPA10mm13cm_P(:,7)==testnum,1)];                     
data13cm_t9 = data13cm_t9(data13cm_t9(:,2)<= maxP &data13cm_t9(:,3)>=minF,:);
data13cm_t10 = data13cm_t10(data13cm_t10(:,2) <= maxP &data13cm_t10(:,3)>=minF,:);
data13cm = [    data13cm_t9(:,1) data13cm_t9(:,2) (data13cm_t9(:,3))/Fmax13; 
                data13cm_t10(:,1) data13cm_t10(:,2) (data13cm_t10(:,3))/Fmax13;
                1 620 0;
                1 620 0];

%23cm
subplot 322
hold on
for a = 1:length(vals_23cm)
    b=vals_23cm(a);
    c = [0 0.424 0.91];
    data23cm_t9 = AllBPA10mm23cm(AllBPA10mm23cm(:,5)==23& AllBPA10mm23cm(:,6)==b&AllBPA10mm23cm(:,7)==testnumber&AllBPA10mm23cm(:,8)==1,:);
    data23cm_t10 = AllBPA10mm23cm(AllBPA10mm23cm(:,5)==23& AllBPA10mm23cm(:,6)==b&AllBPA10mm23cm(:,7)==testnum&AllBPA10mm23cm(:,8)==1,:);
    %txt = sprintf('%dmm',vals_23cm(a));
    txt = sprintf('%d',c(a));
    plot(data23cm_t9(:,2),data23cm_t9(:,1)/Fmax23,'DisplayName',txt)
    plot(data23cm_t10(:,2),data23cm_t10(:,1)/Fmax23,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 23cm all Kinks(Test 9&10)')
data23cm_t9 = [AllBPA10mm23cm_P(AllBPA10mm23cm_P(:,7)==testnumber,13), AllBPA10mm23cm_P(AllBPA10mm23cm_P(:,7)==testnumber,2), AllBPA10mm23cm_P(AllBPA10mm23cm_P(:,7)==testnumber,1)];                     
data23cm_t10 = [AllBPA10mm23cm_P(AllBPA10mm23cm_P(:,7)==testnum,13), AllBPA10mm23cm_P(AllBPA10mm23cm_P(:,7)==testnum,2), AllBPA10mm23cm_P(AllBPA10mm23cm_P(:,7)==testnum,1)];    
data23cm_t9 = data23cm_t9(data23cm_t9(:,2)<= maxP &data23cm_t9(:,3)>=minF,:);
data23cm_t10 = data23cm_t10(data23cm_t10(:,2) <= maxP &data23cm_t10(:,3)>=minF,:);
data23cm = [    data23cm_t9(:,1) data23cm_t9(:,2) (data23cm_t9(:,3))/Fmax23; 
                data23cm_t10(:,1) data23cm_t10(:,2) (data23cm_t10(:,3))/Fmax23;
                1 620 0;
                1 620 0];

%27cm
subplot 323
hold on
for a = 1:length(vals_27cm)
    b = vals_27cm(a);
    c = [0 0.19444 0.41666 0.8611];
    data27cm_t9 = AllBPA10mm27cm(AllBPA10mm27cm(:,5)==27& AllBPA10mm27cm(:,6)==b&AllBPA10mm27cm(:,7)==testnumber&AllBPA10mm27cm(:,8)==1,:);
    data27cm_t10 = AllBPA10mm27cm(AllBPA10mm27cm(:,5)==27& AllBPA10mm27cm(:,6)==b&AllBPA10mm27cm(:,7)==testnum&AllBPA10mm27cm(:,8)==1,:);
    %txt = sprintf('%dmm',vals_27cm(a));
    txt = sprintf('%d',c(a));
    plot(data27cm_t9(:,2),data27cm_t9(:,1)/Fmax27,'DisplayName',txt)
    plot(data27cm_t10(:,2),data27cm_t10(:,1)/Fmax27,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure (kPa)')
ylabel('Force(N)')
title('10mm 27cm all Kinks(Test 9&10)')
data27cm_t9 = [AllBPA10mm27cm_P(AllBPA10mm27cm_P(:,7)==testnumber,13), AllBPA10mm27cm_P(AllBPA10mm27cm_P(:,7)==testnumber,2), AllBPA10mm27cm_P(AllBPA10mm27cm_P(:,7)==testnumber,1)];                     
data27cm_t10 = [AllBPA10mm27cm_P(AllBPA10mm27cm_P(:,7)==testnum,13), AllBPA10mm27cm_P(AllBPA10mm27cm_P(:,7)==testnum,2), AllBPA10mm27cm_P(AllBPA10mm27cm_P(:,7)==testnum,1)];    

data27cm_t9 = data27cm_t9(data27cm_t9(:,2)<= maxP &data27cm_t9(:,3)>=minF,:);
data27cm_t10 = data27cm_t10(data27cm_t10(:,2) <= maxP &data27cm_t10(:,3)>=minF,:);
data27cm = [    data27cm_t9(:,1) data27cm_t9(:,2) (data27cm_t9(:,3))/Fmax27; 
                data27cm_t10(:,1) data27cm_t10(:,2) (data27cm_t10(:,3))/Fmax27;
                1 620 0;
                1 620 0];

%29cm
subplot 324
hold on
for a = 1:length(vals_29cm)
    b = vals_29cm(a);
    c = [0 0.378 0.622 0.911];
    data29cm_t9 = AllBPA10mm29cm(AllBPA10mm29cm(:,5)==29& AllBPA10mm29cm(:,6)==b&AllBPA10mm29cm(:,7)==testnumber&AllBPA10mm29cm(:,8)==1,:);
    data29cm_t10 = AllBPA10mm29cm(AllBPA10mm29cm(:,5)==29& AllBPA10mm29cm(:,6)==b&AllBPA10mm29cm(:,7)==testnum&AllBPA10mm29cm(:,8)==1,:);
    txt = sprintf('%dmm',c(a));
    plot(data29cm_t9(:,2),data29cm_t9(:,1)/Fmax29,'DisplayName',txt)
    plot(data29cm_t10(:,2),data29cm_t10(:,1)/Fmax29,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 29cm all Kinks(Test 9&10)')
data29cm_t9 = [AllBPA10mm29cm_P(AllBPA10mm29cm_P(:,7)==testnumber,13), AllBPA10mm29cm_P(AllBPA10mm29cm_P(:,7)==testnumber,2), AllBPA10mm29cm_P(AllBPA10mm29cm_P(:,7)==testnumber,1)];                     
data29cm_t10 = [AllBPA10mm29cm_P(AllBPA10mm29cm_P(:,7)==testnum,13), AllBPA10mm29cm_P(AllBPA10mm29cm_P(:,7)==testnum,2), AllBPA10mm29cm_P(AllBPA10mm29cm_P(:,7)==testnum,1)];    

data29cm_t9 = data29cm_t9(data29cm_t9(:,2)<= maxP &data29cm_t9(:,3)>=minF,:);
data29cm_t10 = data29cm_t10(data29cm_t10(:,2) <= maxP &data29cm_t10(:,3)>=minF,:);
data29cm = [    data29cm_t9(:,1) data29cm_t9(:,2) (data29cm_t9(:,3))/Fmax29; 
                data29cm_t10(:,1) data29cm_t10(:,2) (data29cm_t10(:,3))/Fmax29;
                1 620 0;
                1 620 0];

%30cm
subplot 325
hold on
for a = 1:length(vals_30cm)
    b = vals_30cm(a);
    c = [0 0.29 0.54 0.80];
    data30cm_t9 = AllBPA10mm30cm(AllBPA10mm30cm(:,5)==30& AllBPA10mm30cm(:,6)==b&AllBPA10mm30cm(:,7)==testnumber&AllBPA10mm30cm(:,8)==1,:);
    data30cm_t10 = AllBPA10mm30cm(AllBPA10mm30cm(:,5)==30& AllBPA10mm30cm(:,6)==b&AllBPA10mm30cm(:,7)==testnum&AllBPA10mm30cm(:,8)==1,:);
    %txt = sprintf('%dmm',vals_30cm(a));
    txt = sprintf('%d',c(a));
    plot(data30cm_t9(:,2),data30cm_t9(:,1)/Fmax30,'DisplayName',txt)
    plot(data30cm_t10(:,2),data30cm_t10(:,1)/Fmax30,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 30cm all Kinks(Test 9&10)')
% for ii = 9:10
%  for j = 1:size(BPA10mm13cm,2)
%    i = ii-8;
%    BPA13cm(j,i) = [BPA10mm13cm{i,j}(BPA10mm13cm{i,j}(:,8)==1,13), (BPA10mm13cm{i,j}(BPA10mm13cm{i,j}(:,8)==1,2)), (BPA10mm13cm{i,j}(BPA10mm13cm{i,j}(:,8)==1,1))];
%  end
% end

data30cm_t9 = [AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnumber,13), AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnumber,2), AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnumber,1)];                     
data30cm_t10 = [AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnum,13), AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnum,2), AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnum,1)];    
data30cm_t9 = data30cm_t9(data30cm_t9(:,2)<= maxP &data30cm_t9(:,3)>=minF,:);
data30cm_t10 = data30cm_t10(data30cm_t10(:,2) <= maxP &data30cm_t10(:,3)>=minF,:);

data30cm = [    data30cm_t9(:,1) data30cm_t9(:,2) (data30cm_t9(:,3))/Fmax30; 
                data30cm_t10(:,1) data30cm_t10(:,2) (data30cm_t10(:,3))/Fmax30;
                1 620 0;
                1 620 0];
            
%% Add Ben's data
%Fmax from experiments
Fmax112 = 325.165;
Fmax415 = 444.8222;
Fmax455 = 479;
Fmax490 = 457.8;
Fmax518 = 429;

%Fmax from maxBPAforce equation
% Fmax112 = maxBPAforce(.112);
% Fmax415 = maxBPAforce(.415);
% Fmax455 = maxBPAforce(.455);
% Fmax490 = maxBPAforce(0.490);
% Fmax518 = maxBPAforce(.518);

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

rawdata45cm = [26	620	0.164835165	1.041666667;
61	620	0.151648352	0.958333333;
150	620	0.101098901	0.638888889;
96	620	0.13956044	0.881944444;
100	620	0.134065934	0.847222222;
119	620	0.12967033	0.819444444;
158	620	0.10989011	0.694444444;
217.0732141	620	0.083516484	0.527777778;
262.0002522	620	0.061538462	0.388888889;
335.3959086	620	0.037362637	0.236111111;
407.0122764	620	0.021978022	0.138888889;
435.9257168	620	0.006593407	0.041666667;
404.7881656	620	0.010989011	0.069444444;
408.3467429	620	0.008791209	0.055555556;
335.3959086	620	0.021978022	0.1388;
257.9968528	620	0.046153846	0.291666667;
203.2837271	620	0.074725275	0.472222222];
data45cm = [rawdata45cm(:,4), rawdata45cm(:,2), rawdata45cm(:,1)/Fmax455];

rawdata49cm = [12	618	0.173469388	0.923913043;
1.11	618	0.175510204	0.934782609;
464	620	-0.002040816	-0.010869565;
433	616	0.008163265	0.043478261;
320	618	0.030612245	0.163043478;
356	620	0.02244898	0.119565217;
300	618	0.03877551	0.206521739;
268	617	0.048979592	0.260869565;
241	617	0.06122449	0.326086957;
215	617	0.069387755	0.369565217;
194	617	0.081632653	0.434782609;
165	617	0.095918367	0.510869565;
144	617	0.106122449	0.565217391;
107	617	0.128571429	0.684782609;
97	617	0.132653061	0.706521739;
75	617	0.146938776	0.782608696;
58	617	0.157142857	0.836956522;
45	617	0.165306122	0.880434783;
20	617	0.17755102	0.945652174;
2.29	617	0.185714286	0.989130435;
0	617	0.187755102	1];
data49cm = [rawdata49cm(:,4), rawdata49cm(:,2), rawdata49cm(:,1)/Fmax490];

rawdata52cm = [0	619	0.166023166	1;
6	618	0.167953668	1.011627907;
-5.3	618	0.171814672	1.034883721;
60	616	0.150579151	0.906976744;
427.6	615	0.005791506	0.034883721;
347	616.8	0.023166023	0.139534884;
296	618	0.038610039	0.23255814;
245	617	0.055984556	0.337209302;
191	617	0.079150579	0.476744186;
55	617	0.152509653	0.918604651;
99	617	0.135135135	0.813953488;
132	617	0.115830116	0.697674419;
160	617	0.102316602	0.61627907;
195	617	0.086872587	0.523255814;
450	617	0.003861004	0.023255814;
420	617	0.011583012	0.069767442;
328	617	0.030888031	0.186046512;
378	617	0.019305019	0.11627907;
230	617	0.063706564	0.38372093;
278	617	0.05019305	0.302325581;
429	617	0	0];
data52cm = [rawdata52cm(:,4), rawdata52cm(:,2), rawdata52cm(:,1)/Fmax518];

%% Combine all data and do a 3d scatter plot. 
allData = [data13cm; data23cm; data27cm; data29cm; data30cm; data11cm; data42cm; data45cm; data49cm; data52cm];
Ben_data = [data11cm; data42cm; data45cm; data49cm; data52cm];

X = allData(:,1); Y=allData(:,2); Z=allData(:,3);
Ynorm = Y/620;

figure
scatter3(X, Y, Z)
xlabel('Relative Strain')
ylabel('Pressure (normalized)')
zlabel('Force (normalized)')
title('10mm force, normalized')

%addnoise to anchor surface at three points
sz = 25;                   % size of anchor surface
err = 0.001;                 % ammount of error
r110 = [1-err + 2*err * rand(sz,2), -err + 2*err * rand(sz,1)].*[1 620 1];
r011 = [-err + 2*err * rand(sz,1), 1-err + 2*err * rand(sz,2)].*[1 620 1];
r000 = -0.01 + 0.02 * rand(sz,3);
anchor = [r110; r011; r000];

%test 9, Big
t9B = [data13cm_t9(:,1) data13cm_t9(:,2) (data13cm_t9(:,3))/Fmax13;
                data23cm_t9(:,1), data23cm_t9(:,2), (data23cm_t9(:,3))/Fmax23; 
                data27cm_t9(:,1) data27cm_t9(:,2) (data27cm_t9(:,3))/Fmax27; 
                data29cm_t9(:,1) data29cm_t9(:,2) (data29cm_t9(:,3))/Fmax29; 
                data30cm_t9(:,1) data30cm_t9(:,2) (data30cm_t9(:,3))/Fmax30]; 
numrows = 400;  
%test 9, resized
t9 = imresize(t9B,[numrows 3]);
                
allData_t9 = [t9;
                Ben_data;
                anchor];

X9 = allData_t9(:,1); Y9=allData_t9(:,2); Z9=allData_t9(:,3);
Y9norm = Y9/620;
            
t10B = [data13cm_t10(:,1) data13cm_t10(:,2) (data13cm_t10(:,3))/Fmax13;
                data23cm_t10(:,1), data23cm_t10(:,2), (data23cm_t10(:,3))/Fmax23; 
                data27cm_t10(:,1) data27cm_t10(:,2) (data27cm_t10(:,3))/Fmax27; 
                data29cm_t10(:,1) data29cm_t10(:,2) (data29cm_t10(:,3))/Fmax29; 
                data30cm_t10(:,1) data30cm_t10(:,2) (data30cm_t10(:,3))/Fmax30]; 
            
t10 = imresize(t10B,[numrows 3]);

allData_t10 = [t10;
                Ben_data;
                anchor];

X10 = t10(:,1); Y10=t10(:,2); Z10=t10(:,3); %Use just test 10 data for validation (i.e. not include anchor or Ben_data)
Y10norm = Y10/620;

figure
hold on
scatter3(X9, Y9norm, Z9, [],'blue')
scatter3(X10, Y10norm, Z10, [],'red')
xlabel('Relative Strain')
ylabel('Pressure (normalized)')
zlabel('Force (normalized)')
title('10mm force, normalized, resized')
hold off

% save allData.mat    allData allData_t9 allData_t10...
%                     t9 t10 Ben_data anchor...
%                     X Y Ynorm Z X9 Y9 Y9norm Z9 X10 Y10 Y10norm Z10 ...
%                     data13cm data23cm data27cm data29cm data30cm data11cm data42cm data45cm data49cm data52cm...
%                     data13cm_t9 data23cm_t9 data27cm_t9 data29cm_t9 data30cm_t9 