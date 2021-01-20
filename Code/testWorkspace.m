theta = deg2rad(linspace(0, 90));
phi = pi/2 - theta/2;

x1 = cos(theta);

for i = 1:length(theta)
    x2(i) = 1 - sin(theta(i))/tan(phi(i));
end

figure
hold on
plot(theta, x1)
plot(theta, x2, '.')
legend('cos', 'sin/tan')
hold off

gamma = linspace(-2*pi, 2*pi);
x3 = tan(gamma)

figure
plot(gamma, x3)
