A = importdata('Joint 1.txt');
dataA = A.data;
time = dataA(:,1);
Ankle_LA = dataA(:,2);
Hip_LA = dataA(:,3);
Knee_LA = dataA(:,4);

B = importdata('Joint 2.txt');
dataB = B.data;
time2 = dataB(:,1);
Ankle_LB = dataB(:,2);
Hip_LB = dataB(:,3);
Knee_LB = dataB(:,4);

Hip_LB(51600:end) = Hip_LA(51600:end);
Knee_LB(56000:end) = Knee_LA(56000:end);
Ankle_LB(70000:end)= Ankle_LA(70000:end);

subplot (3,1,1);
plot(time,Hip_LA,'b',time,Hip_LB,'b:','LineWidth',2);
set(gca,'xticklabel',[],'XMinorTick','on')
ylabel('Hip (Radians)');
% ylim([-0.4 0.3])
grid on
% legend('Ankle','Knee','Hip')

subplot (3,1,2);
plot(time,Knee_LA,'k',time,Knee_LB,'k:','LineWidth',2);
set(gca,'xticklabel',[],'XMinorTick','on')
ylabel('Knee (Radians)');
% ylim([-0.8 0])
grid on
% legend('Ankle','Knee','Hip')

subplot (3,1,3);
plot(time,Ankle_LA,'r',time,Ankle_LB,'r:','LineWidth',2);
set(gca,'XMinorTick','on')
xlabel('Time(s)');
ylabel('Ankle (Radians)');
% ylim([-1 0.5])
grid on
% legend('Ankle','Knee','Hip')
