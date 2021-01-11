p = (-1.2423 + sqrt(1.2423^2 - 4*38.64))/2


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
