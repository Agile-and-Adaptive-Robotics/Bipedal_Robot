
%% Pinned knee, Extensor
%Run and save data from testing results
clear; clc; close all;

Angle = [-102,	-100,	-81.5,	-66,	-47,	-35,	-17,	-11,	-5];
Torque = [5.969270431,	5.781894229,	4.978853364,	3.667219951,	2.623266826,	2.087906249,	1.579313702,	1.311633413,	0];

restingLength = 455;
InflatedLength = [430,	430,	420,	415,	410,	401,	390,	388,	385];
ICRtoMuscle = [30	30	30	30	34	35	45	50	55]/1000;
F = zeros(size(InflatedLength, 2));
TorqueHand = zeros(size(InflatedLength, 2));
for i = 1:size(InflatedLength, 2)
    F(i) = festo2(InflatedLength(i), restingLength, 10, restingLength, 20);
    TorqueHand(i) = ICRtoMuscle(i)*F(i);
end

Theoretical = [13.0706556016290,12.6728468856924,12.2523886598605,11.8096960743697,11.3452250760262,11.0124486010047,10.8962158785483,10.7702424338025,10.6346899078743,10.4897342998847,10.4036056222815,10.3102560360667,10.2095122999010,10.1015106605002,10.0113795417467,9.92367147165927,9.82968609674451,9.72954510571501,9.63506647548198,9.54317145162459,9.44572555393441,9.34284557764291,9.24116839994287,9.14081330643366,9.03546339331611,8.92523502997982,8.81442944799688,8.70395064990682,8.58895347938882,8.46955547636822,8.34827055062158,8.22616103920261,8.09997745642974,7.96984116021790,7.83720248467958,7.70337504163821,7.60015421349621,7.62197526766083,7.63836209591055,7.65113593175467,7.65780953619823,7.65837263135001,7.65346236416486,7.64385866252221,7.62822384770012,7.60656890604511,7.57974446252108,7.54754573834434,7.50944138616414,7.46548757403285,7.41666041255146,7.36207552346272,7.30178437960397,7.23584490192558,7.16432158995732,7.08728566383204,7.00491991138547,6.91755246778993,6.82493738989912,6.72717375117009,6.62391023586357,6.51569852870226,6.40267338328391,6.28456952669727,6.16180368422160,6.03463510260445,5.90267873842748,5.76632293713781,5.62605451028808,5.48159180830614,5.33318086090281,5.18143967479450,5.02583710216024,4.86641815607145,4.70435143206756,4.53955680838773,4.37208043546201,4.25694734679170,4.16769494569301,4.06681213593096,3.95534076443533,3.83304702049298,3.69975653756004,3.55560258718624,3.39879633563622,3.22980352891176,3.04962007733421,2.85769803216360,2.65445259594600,2.43891813979669,2.21056916378079,1.87525877646404,1.71890163926501,1.45552560094059,1.17848926064126,0.888220127493182,0.591624698396069,0.298730802442670,0,0];
phiD = [-119.999999862880,-118.686868548355,-117.373737233831,-116.060605919306,-114.747474604782,-113.434343290257,-112.121211975733,-110.808080661208,-109.494949346683,-108.181818032159,-106.868686717634,-105.555555403110,-104.242424088585,-102.929292774061,-101.616161459536,-100.303030145011,-98.9898988304869,-97.6767675159623,-96.3636362014378,-95.0505048869132,-93.7373735723886,-92.4242422578641,-91.1111109433395,-89.7979796288149,-88.4848483142904,-87.1717169997658,-85.8585856852412,-84.5454543707167,-83.2323230561921,-81.9191917416675,-80.6060604271430,-79.2929291126184,-77.9797977980938,-76.6666664835693,-75.3535351690447,-74.0404038545201,-72.7272725399956,-71.4141412254710,-70.1010099109464,-68.7878785964218,-67.4747472818973,-66.1616159673727,-64.8484846528482,-63.5353533383236,-62.2222220237990,-60.9090907092744,-59.5959593947499,-58.2828280802253,-56.9696967657007,-55.6565654511762,-54.3434341366516,-53.0303028221270,-51.7171715076025,-50.4040401930779,-49.0909088785533,-47.7777775640288,-46.4646462495042,-45.1515149349796,-43.8383836204551,-42.5252523059305,-41.2121209914059,-39.8989896768814,-38.5858583623568,-37.2727270478322,-35.9595957333077,-34.6464644187831,-33.3333331042585,-32.0202017897339,-30.7070704752094,-29.3939391606848,-28.0808078461603,-26.7676765316357,-25.4545452171111,-24.1414139025866,-22.8282825880620,-21.5151512735374,-20.2020199590128,-18.8888886444883,-17.5757573299637,-16.2626260154391,-14.9494947009146,-13.6363633863900,-12.3232320718654,-11.0101007573409,-9.69696944281631,-8.38383812829173,-7.07070681376716,-5.75757549924259,-4.44444418471801,-3.13131287019346,-1.81818155566888,0,0.808081073380252,2.12121238790480,3.43434370242938,4.74747501695393,6.06060633147851,7.37373764600309,8.68686896052764,10.0000002750522];

figure
hold on
title('Isometric Torque vs Knee Angle, Extensor')
xlabel('degrees Flexion(-),Extension(+)')
ylabel('Torque, N*m')
plot(Angle, Torque)
plot(Angle, TorqueHand)
plot(phiD, Theoretical)
legend('FishScale','Moment Arm and Length','Theoretical')
hold off