angle = [10	5	2	0	-5	-10	-15	-30	-45	-60	-75	-90	-105	-120	-135	-150]';
ma = [228.52	228.52	228.52	228.09	227.92	227.76	227.61	227.15	226.63	226	225.28	224.56	223.96	223.65	223.87	225.04]';
fcn = fit(angle, ma, 'cubicspline');

measured_ang = -1*[26.828	31.174	30.936	36.081	39.67	41.189	42.614	45.733	47.827	49.955	53.13	53.673	54.199	56.56	58.392]';
arms = fcn(measured_ang)

Torque_measured = -1*[7.113136537	6.665293601	6.340343593	4.71060257	4.449527719	3.539596694	2.901372937	2.480600054	2.194234049	1.718798005	1.412978615	0.92945505	0.586381564	0.390128872	0.099186634];

figure
plot(phiD, TorqueH(:, 3), phiD, TorqueR(:, 3), measured_ang, Torque_measured, '*')
legend('Human Muscle', 'Calculated BPA', 'Measured BPA')
title('Muscle and PAM Z Torque')
xlabel('Knee Extension/Rotation, degrees')
ylabel('Torque, Nm')
legend('Human', 'BPA Calc', 'BPA measured')
