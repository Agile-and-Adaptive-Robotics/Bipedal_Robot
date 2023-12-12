%% After running LoadAll10mmData and LoadAll10mmData_New from Lawrence's muscle sensory branch
clear allData allData_t9 allData_t10...
                    t9 t10 Ben_data anchor...
                    X Y Ynorm Z X9 Y9 Y9norm Z9 X10 Y10 Y10norm Z10 ...
                    data13cm data23cm data27cm data29cm data30cm...
                    data10cm data15cm data20cm data25cm data30cm_2 data40cm data45cm_2 data52cm_2...
                    data10cm_t9 data15cm_t9 data20cm_t9 data25cm_t9 data30cm_2_t9 data40cm_t9 data45cm_2_t9 data52cm_2_t9...
                    data10cm_t10 data15cm_t10 data20cm_t10 data25cm_t10 data30cm_2_t10 data40cm_t10 data45cm_2_t10 data52cm_2_t10...
                    data11cm data42cm data45cm data49cm data52cm h H...
                    data13cm_t9 data23cm_t9 data27cm_t9 data29cm_t9 data30cm_t9
%% Plot tests 9 and 10, only Pressurizing
 testnumber = 9; %choosing test numner 9 as the standard
 testnum = 10;   %adding test 10 b/c Lawrence says it's good
 

%Mean maximum force at 620kPa, from tests 9 and 10, respecitively.  
Fmax13 = 341.1; %mean([341.9415  340.5389]);
Fmax23 = 377.9; %mean([378.1035 377.6189]);
Fmax27 = 383.4; %mean([383.1724 383.5062]);
Fmax29 = 419.1045; %mean([419.1600 419.0490]);
Fmax30 = 406.9386; %mean([407.2156 406.6616]);

maxP = 800;
minF = 15;

%13cm
figure %11
subplot(5,3,1)
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
subplot(5,3,2)
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
                1 620 0];

%27cm
subplot(5,3,3)
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
                1 620 0];

%29cm
subplot(5,3,4)
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
                1 620 0];

%30cm
subplot(5,3,5)
hold on
for a = 1:length(vals_30cm)
    b = vals_30cm(a);
    c = [0 0.29 0.54 0.80];
    data30cm_t9 = AllBPA10mm30cm(AllBPA10mm30cm(:,5)==30& AllBPA10mm30cm(:,6)==b&AllBPA10mm30cm(:,7)==testnumber&AllBPA10mm30cm(:,8)==1,:);
    data30cm_t10 = AllBPA10mm30cm(AllBPA10mm30cm(:,5)==30& AllBPA10mm30cm(:,6)==b&AllBPA10mm30cm(:,7)==testnum&AllBPA10mm30cm(:,8)==1,:);
    txt = sprintf('%d',c(a));
    plot(data30cm_t9(:,2),data30cm_t9(:,1)/Fmax30,'DisplayName',txt)
    plot(data30cm_t10(:,2),data30cm_t10(:,1)/Fmax30,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 30cm all Kinks(Test 9&10)')

data30cm_t9 = [AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnumber,13), AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnumber,2), AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnumber,1)];                     
data30cm_t10 = [AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnum,13), AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnum,2), AllBPA10mm30cm_P(AllBPA10mm30cm_P(:,7)==testnum,1)];    
data30cm_t9 = data30cm_t9(data30cm_t9(:,2)<= maxP &data30cm_t9(:,3)>=minF,:);
data30cm_t10 = data30cm_t10(data30cm_t10(:,2) <= maxP &data30cm_t10(:,3)>=minF,:);

data30cm = [    data30cm_t9(:,1) data30cm_t9(:,2) (data30cm_t9(:,3))/Fmax30; 
                data30cm_t10(:,1) data30cm_t10(:,2) (data30cm_t10(:,3))/Fmax30;
                1 620 0];
%% Lawrence's 2nd set of accurate tests
Fmax10 = 316.9547325;
Fmax15 = 364;
Fmax20 = 386.2;
Fmax25 = 402.1;
Fmax30_2 = 440.6;
Fmax40 = 465;
Fmax45_2 = 477.7;
Fmax52_2 = 480.6;
% 
% %Make initial guess with maxBPAforce
% % Fmax10 = 300.7;
% % Fmax15 = 361;
% % Fmax20 = 391.5;
% % Fmax25 = 409;
% % Fmax30_2 = 420.5;
% % Fmax40 = 434.6;
% % Fmax45_2 = 438.7042;
% % Fmax52_2 = 445.1;
% 
%10cm
subplot(5,3,6)
hold on
for a = 1:length(kink_p)
    b = kink_p(a);
    c = [8.3,7.975,7.65,7.325];
    data10cm_t9 = AllBPA10mm10cm(AllBPA10mm10cm(:,8)==b &AllBPA10mm10cm(:,14)==testnumber &AllBPA10mm10cm(:,4)==1,:);
    data10cm_t10 = AllBPA10mm10cm(AllBPA10mm10cm(:,8)==b &AllBPA10mm10cm(:,14)==testnum &AllBPA10mm10cm(:,4)==1,:);
    txt = sprintf('%d',b);
    plot(data10cm_t9(:,2),data10cm_t9(:,1)/Fmax10,'DisplayName',txt)
    plot(data10cm_t10(:,2),data10cm_t10(:,1)/Fmax10,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 10cm all Kinks(Test 9&10)')

data10cm_t9 = [AllBPA10mm10cm_P(AllBPA10mm10cm_P(:,14)==testnumber,13), AllBPA10mm10cm_P(AllBPA10mm10cm_P(:,14)==testnumber,2), AllBPA10mm10cm_P(AllBPA10mm10cm_P(:,14)==testnumber,1)];                     
data10cm_t10 = [AllBPA10mm10cm_P(AllBPA10mm10cm_P(:,14)==testnum,13), AllBPA10mm10cm_P(AllBPA10mm10cm_P(:,14)==testnum,2), AllBPA10mm10cm_P(AllBPA10mm10cm_P(:,14)==testnum,1)];    
data10cm_t9 = data10cm_t9(data10cm_t9(:,2)<= maxP &data10cm_t9(:,3)>=minF,:);
data10cm_t10 = data10cm_t10(data10cm_t10(:,2) <= maxP &data10cm_t10(:,3)>=minF,:);
data10cm = [    data10cm_t9(:,1) data10cm_t9(:,2) (data10cm_t9(:,3))/Fmax10; 
                data10cm_t10(:,1) data10cm_t10(:,2) (data10cm_t10(:,3))/Fmax10;
                1 620 0];
            
%15cm
subplot(5,3,7)
hold on
for a = 1:length(kink_p)
    b = kink_p(a);
    c = [13.2	12.675	12.15	11.625];
    data15cm_t9 = AllBPA10mm15cm(AllBPA10mm15cm(:,8)==b&AllBPA10mm15cm(:,14)==testnumber&AllBPA10mm15cm(:,4)==1,:);
    data15cm_t10 = AllBPA10mm15cm(AllBPA10mm15cm(:,8)==b&AllBPA10mm15cm(:,14)==testnum&AllBPA10mm15cm(:,4)==1,:);
    txt = sprintf('%d',c(a));
    plot(data15cm_t9(:,2),data15cm_t9(:,1)/Fmax15,'DisplayName',txt)
    plot(data15cm_t10(:,2),data15cm_t10(:,1)/Fmax15,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 15cm all Kinks(Test 9&10)')

data15cm_t9 = [AllBPA10mm15cm_P((AllBPA10mm15cm_P(:,14)==testnumber),13), AllBPA10mm15cm_P(AllBPA10mm15cm_P(:,14)==testnumber,2), AllBPA10mm15cm_P(AllBPA10mm15cm_P(:,14)==testnumber,1)];                     
data15cm_t10 = [AllBPA10mm15cm_P(AllBPA10mm15cm_P(:,14)==testnum,13), AllBPA10mm15cm_P(AllBPA10mm15cm_P(:,14)==testnum,2), AllBPA10mm15cm_P(AllBPA10mm15cm_P(:,14)==testnum,1)];    
data15cm_t9 = data15cm_t9(data15cm_t9(:,2)<= maxP &data15cm_t9(:,3)>=minF,:);
data15cm_t10 = data15cm_t10(data15cm_t10(:,2) <= maxP &data15cm_t10(:,3)>=minF,:);
data15cm = [    data15cm_t9(:,1) data15cm_t9(:,2) (data15cm_t9(:,3))/Fmax15; 
                data15cm_t10(:,1) data15cm_t10(:,2) (data15cm_t10(:,3))/Fmax15;
                1 620 0];
            
%20cm
subplot(5,3,8)
hold on
for a = 1:length(kink_p)
    b = kink_p(a);
    c = [18.2	17.45	16.7	15.95];
    data20cm_t9 = AllBPA10mm20cm(AllBPA10mm20cm(:,8)==b&AllBPA10mm20cm(:,14)==testnumber&AllBPA10mm20cm(:,4)==1,:);
    data20cm_t10 = AllBPA10mm20cm(AllBPA10mm20cm(:,8)==b&AllBPA10mm20cm(:,14)==testnum&AllBPA10mm20cm(:,4)==1,:);
    txt = sprintf('%d',c(a));
    plot(data20cm_t9(:,2),data20cm_t9(:,1)/Fmax20,'DisplayName',txt)
    plot(data20cm_t10(:,2),data20cm_t10(:,1)/Fmax20,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 20cm all Kinks(Test 9&10)')

data20cm_t9 = [AllBPA10mm20cm_P(AllBPA10mm20cm_P(:,14)==testnumber,13), AllBPA10mm20cm_P(AllBPA10mm20cm_P(:,14)==testnumber,2), AllBPA10mm20cm_P(AllBPA10mm20cm_P(:,14)==testnumber,1)];                     
data20cm_t10 = [AllBPA10mm20cm_P(AllBPA10mm20cm_P(:,14)==testnum,13), AllBPA10mm20cm_P(AllBPA10mm20cm_P(:,14)==testnum,2), AllBPA10mm20cm_P(AllBPA10mm20cm_P(:,14)==testnum,1)];    
data20cm_t9 = data20cm_t9(data20cm_t9(:,2)<= maxP &data20cm_t9(:,3)>=minF,:);
data20cm_t10 = data20cm_t10(data20cm_t10(:,2) <= maxP &data20cm_t10(:,3)>=minF,:);

data20cm = [    data20cm_t9(:,1) data20cm_t9(:,2) (data20cm_t9(:,3))/Fmax20; 
                data20cm_t10(:,1) data20cm_t10(:,2) (data20cm_t10(:,3))/Fmax20;
                1 620 0];
            
%25cm
subplot(5,3,9)
hold on
for a = 1:length(kink_p)
    b = kink_p(a);
    c = [23.3	22.3	21.3	20.3];
    data25cm_t9 = AllBPA10mm25cm(AllBPA10mm25cm(:,8)==b&AllBPA10mm25cm(:,14)==testnumber&AllBPA10mm25cm(:,4)==1,:);
    data25cm_t10 = AllBPA10mm25cm(AllBPA10mm25cm(:,8)==b&AllBPA10mm25cm(:,14)==testnum&AllBPA10mm25cm(:,4)==1,:);
    txt = sprintf('%d',c(a));
    plot(data25cm_t9(:,2),data25cm_t9(:,1)/Fmax25,'DisplayName',txt)
    plot(data25cm_t10(:,2),data25cm_t10(:,1)/Fmax25,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 25cm all Kinks(Test 9&10)')

data25cm_t9 = [AllBPA10mm25cm_P(AllBPA10mm25cm_P(:,14)==testnumber,13), AllBPA10mm25cm_P(AllBPA10mm25cm_P(:,14)==testnumber,2), AllBPA10mm25cm_P(AllBPA10mm25cm_P(:,14)==testnumber,1)];                     
data25cm_t10 = [AllBPA10mm25cm_P(AllBPA10mm25cm_P(:,14)==testnum,13), AllBPA10mm25cm_P(AllBPA10mm25cm_P(:,14)==testnum,2), AllBPA10mm25cm_P(AllBPA10mm25cm_P(:,14)==testnum,1)];    
data25cm_t9 = data25cm_t9(data25cm_t9(:,2)<= maxP &data25cm_t9(:,3)>=minF,:);
data25cm_t10 = data25cm_t10(data25cm_t10(:,2) <= maxP &data25cm_t10(:,3)>=minF,:);

data25cm = [    data25cm_t9(:,1) data25cm_t9(:,2) (data25cm_t9(:,3))/Fmax25; 
                data25cm_t10(:,1) data25cm_t10(:,2) (data25cm_t10(:,3))/Fmax25;
                1 620 0];
            
%30cm_2
subplot(5,3,10)
hold on
for a = 1:length(kink_p)
    b = kink_p(a);
    c = [28.1	26.875	25.65	24.425];
    data30cm_2_t9 = AllBPA10mm30cm_2(AllBPA10mm30cm_2(:,8)==b&AllBPA10mm30cm_2(:,14)==testnumber&AllBPA10mm30cm_2(:,4)==1,:);
    data30cm_2_t10 = AllBPA10mm30cm_2(AllBPA10mm30cm_2(:,8)==b&AllBPA10mm30cm_2(:,14)==testnum&AllBPA10mm30cm_2(:,4)==1,:);
    txt = sprintf('%d',c(a));
    plot(data30cm_2_t9(:,2),data30cm_2_t9(:,1)/Fmax30_2,'DisplayName',txt)
    plot(data30cm_2_t10(:,2),data30cm_2_t10(:,1)/Fmax30_2,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 30cm_2 all Kinks(Test 9&10)')

data30cm_2_t9 = [AllBPA10mm30cm_2_P(AllBPA10mm30cm_2_P(:,14)==testnumber,13), AllBPA10mm30cm_2_P(AllBPA10mm30cm_2_P(:,14)==testnumber,2), AllBPA10mm30cm_2_P(AllBPA10mm30cm_2_P(:,14)==testnumber,1)];                     
data30cm_2_t10 = [AllBPA10mm30cm_2_P(AllBPA10mm30cm_2_P(:,14)==testnum,13), AllBPA10mm30cm_2_P(AllBPA10mm30cm_2_P(:,14)==testnum,2), AllBPA10mm30cm_2_P(AllBPA10mm30cm_2_P(:,14)==testnum,1)];    
data30cm_2_t9 = data30cm_2_t9(data30cm_2_t9(:,2)<= maxP &data30cm_2_t9(:,3)>=minF,:);
data30cm_2_t10 = data30cm_2_t10(data30cm_2_t10(:,2) <= maxP &data30cm_2_t10(:,3)>=minF,:);

data30cm_2 = [ data30cm_2_t9(:,1) data30cm_2_t9(:,2) (data30cm_2_t9(:,3))/Fmax30_2; 
                data30cm_2_t10(:,1) data30cm_2_t10(:,2) (data30cm_2_t10(:,3))/Fmax30_2;
                1 620 0];

%40cm
subplot(5,3,11)
hold on
for a = 1:length(kink_p)
    b = kink_p(a);
    c = [38.2	36.5	34.8	33.1];
    data40cm_t9 = AllBPA10mm40cm(AllBPA10mm40cm(:,8)==b&AllBPA10mm40cm(:,14)==testnumber&AllBPA10mm40cm(:,4)==1,:);
    data40cm_t10 = AllBPA10mm40cm(AllBPA10mm40cm(:,8)==b&AllBPA10mm40cm(:,14)==testnum&AllBPA10mm40cm(:,4)==1,:);
    txt = sprintf('%d',c(a));
    plot(data40cm_t9(:,2),data40cm_t9(:,1)/Fmax40,'DisplayName',txt)
    plot(data40cm_t10(:,2),data40cm_t10(:,1)/Fmax40,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 40cm all Kinks(Test 9&10)')

data40cm_t9 = [AllBPA10mm40cm_P(AllBPA10mm40cm_P(:,14)==testnumber,13), AllBPA10mm40cm_P(AllBPA10mm40cm_P(:,14)==testnumber,2), AllBPA10mm40cm_P(AllBPA10mm40cm_P(:,14)==testnumber,1)];                     
data40cm_t10 = [AllBPA10mm40cm_P(AllBPA10mm40cm_P(:,14)==testnum,13), AllBPA10mm40cm_P(AllBPA10mm40cm_P(:,14)==testnum,2), AllBPA10mm40cm_P(AllBPA10mm40cm_P(:,14)==testnum,1)];    
data40cm_t9 = data40cm_t9(data40cm_t9(:,2)<= maxP &data40cm_t9(:,3)>=minF,:);
data40cm_t10 = data40cm_t10(data40cm_t10(:,2) <= maxP &data40cm_t10(:,3)>=minF,:);

data40cm = [ data40cm_t9(:,1) data40cm_t9(:,2) (data40cm_t9(:,3))/Fmax40; 
                data40cm_t10(:,1) data40cm_t10(:,2) (data40cm_t10(:,3))/Fmax40;
                1 620 0];
            
%45cm, use _2 to not confuse with other 40cm, unkinked only
subplot(5,3,12)
hold on
for a = 1
    b = kink_p(a);
    c = [42.6];
    data45cm_2_t9 = AllBPA10mm45cm_2(AllBPA10mm45cm_2(:,8)==b&AllBPA10mm45cm_2(:,14)==testnumber&AllBPA10mm45cm_2(:,4)==1,:);
    data45cm_2_t10 = AllBPA10mm45cm_2(AllBPA10mm45cm_2(:,8)==b&AllBPA10mm45cm_2(:,14)==testnum&AllBPA10mm45cm_2(:,4)==1,:);
    txt = sprintf('%d',c(a));
    plot(data45cm_2_t9(:,2),data45cm_2_t9(:,1)/Fmax45_2,'DisplayName',txt)
    plot(data45cm_2_t10(:,2),data45cm_2_t10(:,1)/Fmax45_2,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 45cm_2 all Kinks(Test 9&10)')

data45cm_2_t9 = [AllBPA10mm45cm_2_P(AllBPA10mm45cm_2_P(:,14)==testnumber,13), AllBPA10mm45cm_2_P(AllBPA10mm45cm_2_P(:,14)==testnumber,2), AllBPA10mm45cm_2_P(AllBPA10mm45cm_2_P(:,14)==testnumber,1)];                     
data45cm_2_t10 = [AllBPA10mm45cm_2_P(AllBPA10mm45cm_2_P(:,14)==testnum,13), AllBPA10mm45cm_2_P(AllBPA10mm45cm_2_P(:,14)==testnum,2), AllBPA10mm45cm_2_P(AllBPA10mm45cm_2_P(:,14)==testnum,1)];    
data45cm_2_t9 = data45cm_2_t9(data45cm_2_t9(:,2)<= maxP &data45cm_2_t9(:,3)>=minF,:);
data45cm_2_t10 = data45cm_2_t10(data45cm_2_t10(:,2) <= maxP &data45cm_2_t10(:,3)>=minF,:);

data45cm_2 = [ data45cm_2_t9(:,1) data45cm_2_t9(:,2) (data45cm_2_t9(:,3))/Fmax45_2; 
                data45cm_2_t10(:,1) data45cm_2_t10(:,2) (data45cm_2_t10(:,3))/Fmax45_2;
                1 620 0];
            
%52cm_2, use _2 to not confuse with other 52cm BPA
subplot(5,3,13)
hold on
for a = 1:length(kink_p)
    b = kink_p(a);
    c = [52.1 50 47.9 45.8];
    data52cm_2_t9 = AllBPA10mm52cm_2_P(AllBPA10mm52cm_2_P(:,8)==b&AllBPA10mm52cm_2_P(:,14)==testnumber&AllBPA10mm52cm_2_P(:,4)==1,:);
    data52cm_2_t10 = AllBPA10mm52cm_2_P(AllBPA10mm52cm_2_P(:,8)==b&AllBPA10mm52cm_2_P(:,14)==testnum&AllBPA10mm52cm_2_P(:,4)==1,:);
    txt = sprintf('%d',c(a));
    plot(data52cm_2_t9(:,2),data52cm_2_t9(:,1)/Fmax52_2,'DisplayName',txt)
    plot(data52cm_2_t10(:,2),data52cm_2_t10(:,1)/Fmax52_2,'DisplayName',txt)
end
hold off
legend
xlabel('Pressure(kPa)')
ylabel('Force(N)')
title('10mm 52cm_2 all Kinks(Test 9&10)')

data52cm_2_t9 = [AllBPA10mm52cm_2_P(AllBPA10mm52cm_2_P(:,14)==testnumber,13), AllBPA10mm52cm_2_P(AllBPA10mm52cm_2_P(:,14)==testnumber,2), AllBPA10mm52cm_2_P(AllBPA10mm52cm_2_P(:,14)==testnumber,1)];                     
data52cm_2_t10 = [AllBPA10mm52cm_2_P(AllBPA10mm52cm_2_P(:,14)==testnum,13), AllBPA10mm52cm_2_P(AllBPA10mm52cm_2_P(:,14)==testnum,2), AllBPA10mm52cm_2_P(AllBPA10mm52cm_2_P(:,14)==testnum,1)];    
data52cm_2_t9 = data52cm_2_t9(data52cm_2_t9(:,2)<= maxP &data52cm_2_t9(:,3)>=minF,:);
data52cm_2_t10 = data52cm_2_t10(data52cm_2_t10(:,2) <= maxP &data52cm_2_t10(:,3)>=minF,:);

data52cm_2 = [ data52cm_2_t9(:,1) data52cm_2_t9(:,2) (data52cm_2_t9(:,3))/Fmax52_2; 
                data52cm_2_t10(:,1) data52cm_2_t10(:,2) (data52cm_2_t10(:,3))/Fmax52_2;
                1 620 0];
%% Add Ben's data
%Fmax from experiments, using interpolation
Fmax112 = 334.74;
Fmax415 = 445;
Fmax455 = 460;
Fmax490 = 458.6;
Fmax518 = 455.83;

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
allData = [data13cm; data23cm; data27cm; data29cm; data30cm; ...
            data10cm; data15cm; data20cm; data25cm; data30cm_2; data40cm; data45cm_2; data52cm_2;...
            data11cm; data42cm; data45cm; data49cm; data52cm];
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
sz = 15;                   % size of anchor surface
err = 0.005;                 % ammount of error
r110 = [1-err + 2*err * rand(sz,2), -err + 2*err * rand(sz,1)].*[1 620 1];
r011 = [-err + 2*err * rand(sz,1), 1-err + 2*err * rand(sz,2)].*[1 620 1];
r000 = -0.01 + 0.02 * rand(sz,3);
anchor = [r110; r011; r000];

%test 9, Big
t9B = [data13cm_t9(:,1) data13cm_t9(:,2) (data13cm_t9(:,3))/Fmax13;
                data23cm_t9(:,1), data23cm_t9(:,2), (data23cm_t9(:,3))/Fmax23; 
                data27cm_t9(:,1) data27cm_t9(:,2) (data27cm_t9(:,3))/Fmax27; 
                data29cm_t9(:,1) data29cm_t9(:,2) (data29cm_t9(:,3))/Fmax29; 
                data30cm_t9(:,1) data30cm_t9(:,2) (data30cm_t9(:,3))/Fmax30;
%                 data10cm_t9(:,1) data10cm_t9(:,2) (data10cm_t9(:,3))/Fmax10;
%                 data15cm_t9(:,1) data15cm_t9(:,2) (data15cm_t9(:,3))/Fmax15;
%                 data20cm_t9(:,1) data20cm_t9(:,2) (data20cm_t9(:,3))/Fmax20;
%                 data25cm_t9(:,1) data25cm_t9(:,2) (data25cm_t9(:,3))/Fmax25;
%                 data30cm_2_t9(:,1) data30cm_2_t9(:,2) (data30cm_2_t9(:,3))/Fmax30_2;
%                 data40cm_t9(:,1) data40cm_t9(:,2) (data40cm_t9(:,3))/Fmax40;
%                 data45cm_2_t9(:,1) data45cm_2_t9(:,2) (data45cm_2_t9(:,3))/Fmax45_2;
%                 data52cm_2_t9(:,1) data52cm_2_t9(:,2) (data52cm_2_t9(:,3))/Fmax52_2;
                ]; 
numrows = 200;  
%test 9, resized
law = {data13cm_t9; data23cm_t9; data27cm_t9; data29cm_t9; data30cm_t9};
sm = numrows/length(law);

for i=1:length(law)
   h{i,1} = imresize(law{i},[sm 3],'nearest');
end
t9 = cell2mat(h);

allData_t9 = [t9B;
                Ben_data;
                anchor];

X9 = allData_t9(:,1); Y9=allData_t9(:,2); Z9=allData_t9(:,3);
Y9norm = Y9/620;
            
t10B = [data13cm_t10(:,1) data13cm_t10(:,2) (data13cm_t10(:,3))/Fmax13;
                data23cm_t10(:,1), data23cm_t10(:,2), (data23cm_t10(:,3))/Fmax23; 
                data27cm_t10(:,1) data27cm_t10(:,2) (data27cm_t10(:,3))/Fmax27; 
                data29cm_t10(:,1) data29cm_t10(:,2) (data29cm_t10(:,3))/Fmax29; 
                data30cm_t10(:,1) data30cm_t10(:,2) (data30cm_t10(:,3))/Fmax30;
%                 data10cm_t10(:,1) data10cm_t10(:,2) (data10cm_t10(:,3))/Fmax10;
%                 data15cm_t10(:,1) data15cm_t10(:,2) (data15cm_t10(:,3))/Fmax15;
%                 data20cm_t10(:,1) data20cm_t10(:,2) (data20cm_t10(:,3))/Fmax20;
%                 data25cm_t10(:,1) data25cm_t10(:,2) (data25cm_t10(:,3))/Fmax25;
%                 data30cm_2_t10(:,1) data30cm_2_t10(:,2) (data30cm_2_t10(:,3))/Fmax30_2;
%                 data40cm_t10(:,1) data40cm_t10(:,2) (data40cm_t10(:,3))/Fmax40;
%                 data45cm_2_t10(:,1) data45cm_2_t10(:,2) (data45cm_2_t10(:,3))/Fmax45_2;
%                 data52cm_2_t10(:,1) data52cm_2_t10(:,2) (data52cm_2_t10(:,3))/Fmax52_2;
                ]; 
            
%test 10, resized

LAW = {data13cm_t10; data23cm_t10; data27cm_t10; data29cm_t10; data30cm_t10};
sm = numrows/length(LAW);

for i=1:length(law)
   H{i,1} = imresize(LAW{i},[sm 3],'nearest');
end
t10 = cell2mat(H);

allData_t10 = [t10B;
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

save allData.mat    allData allData_t9 allData_t10...
                    t9 t10 Ben_data anchor...
                    X Y Ynorm Z X9 Y9 Y9norm Z9 X10 Y10 Y10norm Z10 ...
                    data13cm data23cm data27cm data29cm data30cm...
                    data11cm data42cm data45cm data49cm data52cm...
                    data13cm_t9 data23cm_t9 data27cm_t9 data29cm_t9 data30cm_t9 ...
                    data10cm_t9 data15cm_t9 data20cm_t9 data25cm_t9 data30cm_2_t9 data40cm_t9 data45cm_2_t9 data52cm_2_t9...