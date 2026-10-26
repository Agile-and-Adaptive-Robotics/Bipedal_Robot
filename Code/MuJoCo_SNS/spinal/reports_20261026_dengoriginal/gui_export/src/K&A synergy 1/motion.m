A = importdata('Joint.txt');
data = A.data;
time = data(:,1);
Ankle_L = data(:,2);
Hip_L = data(:,3);
Knee_L = data(:,5);

subplot (3,1,1);
plot(time,Hip_L,'b','LineWidth',2);
set(gca,'xticklabel',[],'XMinorTick','on')
ylabel('Hip (Radians)');
% ylim([-0.4 0.3])
grid on
% legend('Ankle','Knee','Hip')

subplot (3,1,2);
plot(time,Knee_L,'k','LineWidth',2);
set(gca,'xticklabel',[],'XMinorTick','on')
ylabel('Knee (Radians)');
% ylim([-0.8 0])
grid on
% legend('Ankle','Knee','Hip')

subplot (3,1,3);
plot(time,Ankle_L,'r','LineWidth',2);
set(gca,'XMinorTick','on')
xlabel('Time(s)');
ylabel('Ankle (Radians)');
% ylim([-1 0.5])
grid on
% legend('Ankle','Knee','Hip')
