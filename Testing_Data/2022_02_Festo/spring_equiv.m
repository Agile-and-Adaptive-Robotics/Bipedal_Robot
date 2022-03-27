k_BPA = 509; 
kPSE = linspace(0, 100);

kEQ_par = zeros(1,size(kPSE,2));
kEQ_ser = zeros(1,size(kPSE,2));
for i = 1:size(kPSE,2)
    kEQ_par(i) = (k_BPA*kPSE(i))/(k_BPA+kPSE(i));
    kEQ_ser(i) = k_BPA+kPSE(i);
end

figure
plot(kPSE, kEQ_par, kPSE, kEQ_ser)
xlabel('Spring Constant of Parallel or Series Element')
ylabel('Spring Equivalent of System')
title('Affect of added spring to system')
legend('Parallel','Series')
