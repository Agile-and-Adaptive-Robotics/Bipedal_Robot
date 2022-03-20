%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

restingLength = 0.485; %resting length, m
kmax = 0.405; %Length at maximum contraction, m

load KneeExtPin_10mm_48cm.mat
Theoretical = TorqueR(:,3)';
phiD = [-119.999999862880,-118.686868548355,-117.373737233831,-116.060605919306,-114.747474604782,-113.434343290257,-112.121211975733,-110.808080661208,-109.494949346683,-108.181818032159,-106.868686717634,-105.555555403110,-104.242424088585,-102.929292774061,-101.616161459536,-100.303030145011,-98.9898988304869,-97.6767675159623,-96.3636362014378,-95.0505048869132,-93.7373735723886,-92.4242422578641,-91.1111109433395,-89.7979796288149,-88.4848483142904,-87.1717169997658,-85.8585856852412,-84.5454543707167,-83.2323230561921,-81.9191917416675,-80.6060604271430,-79.2929291126184,-77.9797977980938,-76.6666664835693,-75.3535351690447,-74.0404038545201,-72.7272725399956,-71.4141412254710,-70.1010099109464,-68.7878785964218,-67.4747472818973,-66.1616159673727,-64.8484846528482,-63.5353533383236,-62.2222220237990,-60.9090907092744,-59.5959593947499,-58.2828280802253,-56.9696967657007,-55.6565654511762,-54.3434341366516,-53.0303028221270,-51.7171715076025,-50.4040401930779,-49.0909088785533,-47.7777775640288,-46.4646462495042,-45.1515149349796,-43.8383836204551,-42.5252523059305,-41.2121209914059,-39.8989896768814,-38.5858583623568,-37.2727270478322,-35.9595957333077,-34.6464644187831,-33.3333331042585,-32.0202017897339,-30.7070704752094,-29.3939391606848,-28.0808078461603,-26.7676765316357,-25.4545452171111,-24.1414139025866,-22.8282825880620,-21.5151512735374,-20.2020199590128,-18.8888886444883,-17.5757573299637,-16.2626260154391,-14.9494947009146,-13.6363633863900,-12.3232320718654,-11.0101007573409,-9.69696944281631,-8.38383812829173,-7.07070681376716,-5.75757549924259,-4.44444418471801,-3.13131287019346,-1.81818155566888,0,0.808081073380252,2.12121238790480,3.43434370242938,4.74747501695393,6.06060633147851,7.37373764600309,8.68686896052764,10.0000002750522];
%% Tests 1 and 4 done with CALT load cell. Tests 2 and 3 done with fish scale. Fish scale tests had pressure spot checked around 612 kPa. 


%% Torque calculated from measurements

Angle1 = [-121.5	-108.5	-101.5	-94	-83	-75	-61.5	-65	-32	-18	-32.5	-47.5	-57.5	-69.5	-80	-90.5	-96.5	-104	-113.5	-123];
Angle2 = [-110	-82	-71.5	-55	-42	-43	-65.5	-74.5	-72	-82	-88.5	-98	-113.5];
Angle3 = [-110	-105	-103.5	-91	-80	-72.5	-62	-48	-40.5	-31	-30	-31.5	-43	-50	-62.5	-69	-78	-88	-93	-100	-106	-107];
Angle4 = [-123	-119	-106	-92];
Angle = [Angle1, Angle2, Angle3, Angle4];

Torque1 = [4.62685036	4.503204046	4.29202985	4.002027833	3.417014674	2.826201474	1.538232508	1.214084497	0.439203694	-0.080856656	0.539217698	1.260349872	2.254964352	3.207688111	3.851548806	4.243514283	4.672306017	4.873564883	5.146026298	6.087161098];
Torque2 = [5.783034531	4.71210221	3.667943198	1.847358253	1.044159012	1.579625173	2.356051105	3.159250345	3.50730335	4.203409358	4.203409358	4.97983529	4.97983529];
Torque3 = [5.783034531	4.71210221	4.471142438	3.935676278	3.667943198	3.159250345	2.623784185	1.847358253	1.419266981	0.53546616	0	0.53546616	1.044159012	2.356051105	2.891517265	3.400210118	3.935676278	4.203409358	4.71210221	4.97983529	6.291727383	6.559460463];
Torque4 = [7.013441569	4.72935096	4.346134571	3.454623635];
Torque = [Torque1, Torque2, Torque3, Torque4];

%% Calculate Torque by finding force from muscle contraction and distance
%from force line of action to muscle ICR.
%Hand measurements for Fish scale Extensor test 2 were done incorrectly and will be disregarded

InflatedLength1 = [455.5	461	451	451	440	438.5	435	415.5	413	410	409.5	418	426.5	438	442.5	441.5	450	458.5	455.5	463]/1000;
%InflatedLength2 = [448	445.5	429.5	413.5	401	410	417.5	422.5	430.5	429.5	441.5	439	443.5]/1000;
InflatedLength3 = [445	445	433	426	423	421	420	415	405	400	398	400	405	407	416	420	430	433	441	440	450	453]/1000;
InflatedLength4 = [459	456	447	438]/1000;
InflatedLength = [InflatedLength1, InflatedLength3, InflatedLength4];

ICRtoMuscle1 = [30.5	41.5	39	35.5	37	40.5	37	46.5	38.5	45.5	36.5	38.5	35	38.5	36	42	39	36	38	37]/1000;
%ICRtoMuscle2 = [46	54.5	53	52	57	49	53	53	54	51	52	53	55.5]/1000;
ICRtoMuscle3 = [32	30	30	30	30	30	30	30	35	35	40	35	34	33	30	30	30	30	30	30	30	30]/1000;
ICRtoMuscle4 = [30	30	30	30]/1000;
ICRtoMuscle = [ICRtoMuscle1, ICRtoMuscle3, ICRtoMuscle4];

F1 = zeros(1,size(InflatedLength1, 2));
%F2 = zeros(1,size(InflatedLength2, 2));
F3 = zeros(1,size(InflatedLength3, 2));
F4 = zeros(1,size(InflatedLength4, 2));
F = [F1, F3, F4];

TorqueHand1 = zeros(1,size(InflatedLength1, 2));
%TorqueHand2 = zeros(1,size(InflatedLength2, 2));
TorqueHand3 = zeros(1,size(InflatedLength3, 2));
TorqueHand4 = zeros(1,size(InflatedLength4, 2));
TorqueHand = [TorqueHand1, TorqueHand3, TorqueHand4];

%load pressure where applicable
test = [1 3];
runsperseries = [20 4];

    pres1 = zeros(1,size(runsperseries(1),2));
    pres3 = zeros(1,size(runsperseries(2),2));
    
    for j = 1:runsperseries(1)
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(1),j)
                load(file_name,'Stats')
                pres1(1,j) = Stats{'Mean',2}
    end
    for j = 1:runsperseries(2)
                file_name = sprintf('ExtTest%0.0f_%0.0f.mat', test(2),j)
                load(file_name,'Stats')
                pres3(1,j) = Stats{'Mean',2}
    end
pres2 = 612*ones(1,size(InflatedLength3, 2));
pres = [pres1 pres2 pres3];

for i = 1:size(InflatedLength, 2)
    F(i) = festo3(InflatedLength(i), restingLength, 10, pres(i), 0.405);    
    TorqueHand(i) = ICRtoMuscle(i)*F(i);
end
TorqueHand1 = TorqueHand(1:size(TorqueHand1,2));
TorqueHand3 = TorqueHand(((size(TorqueHand1,2)+1)):(size(TorqueHand1,2)+size(TorqueHand3,2)));
TorqueHand4 = TorqueHand((((size(TorqueHand1,2)+size(TorqueHand3,2))+1)):size(TorqueHand,2));
%% Mean and RMSE
X1 = linspace(min(Angle),max(Angle),size(Angle,2));      %Range of motion
X2 = linspace(min(Angle),max(Angle),(size(Angle1,2)+size(Angle3,2)+size(Angle4,2)));

modp = 'poly3';
fitOp = fitoptions(modp,'Normalize','on');
[mdl1, gofp1] = fit(Angle',Torque',modp,fitOp)
TorqueStd = gofp1.rmse
TorqueMean = feval(mdl1,X1)';

[mdl2, gofp2] = fit([Angle1 Angle3 Angle4]',TorqueHand',modp,fitOp);
HandStd = gofp2.rmse;
HandMean = feval(mdl2,X2)';

mod = fittype( {'(sind(x-1))', '((x-10)^2)', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b', 'c'} );
fitOptions = fitoptions(mod);
[mdl1u, gof1] = fit(Angle',Torque',mod,fitOptions)
TorqueStdu = gof1.rmse
TorqueMeanu = feval(mdl1u,X1)';

[mdl2u, gof2] = fit([Angle1 Angle3 Angle4]',TorqueHand',mod,fitOptions);
HandStdu = gof2.rmse;
HandMeanu = feval(mdl2u,X2)';

%% Plotting the polynomial solution
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Extensor, 48.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
plot(phiD, Theoretical,'Color',[0 0.4470 0.7410],'Linewidth',2,'DisplayName','Theoretical Calculation')

X1new=[X1,fliplr(X1)];
Y1=[TorqueMean+TorqueStd,fliplr(TorqueMean-TorqueStd)];
X2new=[X2,fliplr(X2)];
Y2=[HandMean+HandStd,fliplr(HandMean-HandStd)];
plot(X1,TorqueMean,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
fill(X1new,Y1,[1 0.4 0.8],'DisplayName','Scale SD','FaceAlpha',0.25);
plot(X2,HandMean,'--r','Linewidth',2,'DisplayName','Torque mean, hand')
fill(X2new,Y2,[.6 1.0 .6],'DisplayName','Hand torque SD','FaceAlpha',0.25);

sz = 50;
c = [0.8500 0.3250 0.0980]; % color
scatter(Angle1,Torque1,sz,'d','g','DisplayName','BB LC');
scatter(Angle2,Torque2,sz,'d','r','DisplayName','JM FS');
scatter(Angle3,Torque3,sz,'d','b','DisplayName','BB FS');
scatter(Angle4,Torque4,sz,'d','CData',c,'DisplayName','BB LC');
scatter(Angle1,TorqueHand1,[],'g','filled','DisplayName','JM hand');
%scatter(Angle2,TorqueHand2,[],'r','filled','DisplayName','JM hand');
scatter(Angle3,TorqueHand3,[],'b','filled','DisplayName','BB hand');
scatter(Angle4,TorqueHand4,[],[0.8500 0.3250 0.0980],'DisplayName','BB hand');

legend
hold off
%% Plotting nonlinear solution
figure
hold on
title('Isometric Torque vs Knee Angle, 10mm Extensor, 48.5cm long')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
plot(phiD, Theoretical,'Color',[0 0.4470 0.7410],'Linewidth',2,'DisplayName','Theoretical Calculation')

Y3=[TorqueMeanu+TorqueStdu,fliplr(TorqueMeanu-TorqueStdu)];
Y4=[HandMeanu+HandStdu,fliplr(HandMeanu-HandStdu)];
plot(X1,TorqueMeanu,'--k','Linewidth',2,'DisplayName','Torque mean, scale')
fill(X1new,Y3,[1 0.4 0.8],'DisplayName','Scale SD','FaceAlpha',0.25);
plot(X2,HandMeanu,'--r','Linewidth',2,'DisplayName','Torque mean, hand')
fill(X2new,Y4,[.6 1.0 .6],'DisplayName','Hand torque SD','FaceAlpha',0.25);

sz = 50;
c = [0.8500 0.3250 0.0980]; % color
scatter(Angle1,Torque1,sz,'d','g','DisplayName','BB LC');
scatter(Angle2,Torque2,sz,'d','r','DisplayName','JM FS');
scatter(Angle3,Torque3,sz,'d','b','DisplayName','BB FS');
scatter(Angle4,Torque4,sz,'d','CData',c,'DisplayName','BB LC');
scatter(Angle1,TorqueHand1,[],'g','filled','DisplayName','JM hand');
%scatter(Angle2,TorqueHand2,[],'r','filled','DisplayName','JM hand');
scatter(Angle3,TorqueHand3,[],'b','filled','DisplayName','BB hand');
scatter(Angle4,TorqueHand4,[],[0.8500 0.3250 0.0980],'DisplayName','BB hand');

legend
hold off