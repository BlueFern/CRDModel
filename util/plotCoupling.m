% Plot the coupling strength vs theta

clear;

theta = linspace(0,2*pi,500);
R1 = 80/(2*pi);
R2 = 40/(2*pi);
r = 20/(2*pi);
a1 = sqrt(R1^2 - r^2)
a2 = sqrt(R2^2 - r^2)
eta1 = atanh(a1/R1)
eta2 = atanh(a2/R2)

theta_i1 = acos(R1/r - (a1^2)./(r*(R1+r*cos(theta))));
theta_i2 = acos(R2/r - (a2^2)./(r*(R2+r*cos(theta))));

C1 = 10*(cosh(eta1)-cos(theta_i1)).^2/(a1^2);
C2 = 10*(cosh(eta2)-cos(theta_i2)).^2/(a2^2);

%hold on
figure(1);
plot(theta, C2, theta, C1, [0, 2*pi], [1,1]);
title('Coupling Strength');
xlim([0 2*pi]);
%hold off
legend('R=40','R=80','Flat');
xlabel('theta')