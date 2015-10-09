% Plot the coupling strength vs theta

clear;

theta = linspace(0,pi,100);
R = 40/(2*pi)
r = 20/(2*pi)
a = sqrt(R^2 - r^2)
eta = atanh(a/R)

theta_i = acos(R/r - (a^2)./(r*(R+r*cos(theta))));

C = (cosh(eta)-cos(theta_i)).^2/(a^2);

hold on
plot(theta, C);
title('Coupling Strength')
hold off