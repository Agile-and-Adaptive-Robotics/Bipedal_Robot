sigM1 = (10 - 20)/2

sigA1 = (10 + 20)/2

sigRev1 = sigA1/(1 - sigM1/200)

a = (0.77 * 200)^2/25

b = -1/3*log10(0.77*200/25)

N1 = (sigRev1/a)^(1/b)

sigM2 = (40 - 20)/2

sigA2 = (40 + 20)/2

sigRev2 = sigA2/(1 - sigM2/200)

N2 = (sigRev2/a)^(1/b)

sigM3 = (60 - 40)/2

sigA3 = (60 + 40)/2

sigRev3 = sigA3/(1 - sigM3/200)

N3 = (sigRev3/a)^(1/b)

D = 2/N1 + 3/N2 + 2/N3

t = 6/D/60/60


% dH = 0.25/1.5
% 
% sigMax = 12/(0.25*1.25)
% 
% sigMax2 = 12000/(0.25)
% 
% sigMin = -5/(0.25*1.25)
% 
% sigMin2 = -5000/0.25
% 
% sigM1 = (sigMax + sigMin)/2
% 
% sigM2 = (sigMax2 + sigMin2)/2
% 
% sigA1 = (sigMax - sigMin)/2
% 
% sigA2 = (sigMax2 - sigMin2)/2
% 
% Kfh = 1 + 0.95*(2.6-1)
% 
% Kff = 1 + 0.95*(1.8 - 1)
% 
% sigM1 = Kfh*sigM1
% 
% sigA1 = Kfh*sigA1
% 
% sigM2 = Kff*sigM2
% 
% sigA2 = Kff*sigA2
% 
% sigRev = sigA1/(1 - sigM1/185)
% 
% a = (0.78 * 185)^2/40
% 
% b = -1/3*log10(0.78*185/40)
% 
% N = (sigRev/a)^(1/b)



% clc                 %Command Line Clear
% clear               %Clear the workspace of stored variables
% close all           %Close all open figures
% 
% x = linspace(0, 2);
% 
% muscle = 3*x - 3;
% 
% Pam1 = 4*x + 7;
% Pam1 = 8*ones(100, 1);
% 
% Pam2 = 4*x - 4;
% 
% figure
% hold on
% plot(x, muscle)
% plot(x, Pam1)
% plot(x, Pam2)
% plot(x, linspace(0, 0), 'k')
% legend('Muscle', 'Pam1', 'Pam2')
% hold off
% 
% Pam1C = 0;
% Pam2C = 0;
% for i = 1:length(x)
%     if muscle(i) >= 0
%         Pam1C = Pam1C + Pam1(i) - muscle(i);
%         Pam2C = Pam2C + Pam2(i) - muscle(i);
%     elseif muscle(i) < 0
%         Pam1C = Pam1C + muscle(i) - Pam1(i);
%         Pam2C = Pam2C + muscle(i) - Pam2(i);
%     end
% end
% 
% 
% muscleAbs = abs(muscle);
% Pam2Abs = abs(Pam2);
% 
% figure
% hold on
% plot(x, muscleAbs)
% plot(x, Pam1)
% plot(x, Pam2Abs)
% plot(x, linspace(0, 0), 'k')
% legend('Muscle', 'Pam1', 'Pam2')
% hold off
% 
% muscleNom = muscle./muscle;
% Pam1Nom = Pam1./muscle;
% Pam2Nom = Pam2./muscle;
% 
% figure
% hold on
% plot(x, muscleNom)
% plot(x, Pam1Nom)
% plot(x, Pam2Nom)
% plot(x, linspace(0, 0), 'k')
% legend('Muscle', 'Pam1', 'Pam2')
% hold off
% 
% siga = 1.63*2.5
% sigm = 1.63*12.5
