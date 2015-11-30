% Plot the Gaussian Curvature vs theta

clear;

theta = linspace(0,2*pi,500);
R1 = 80/(2*pi);
R2 = 40/(2*pi);
r = 20/(2*pi);

G1 = cos(theta)./(r*(R1+r*cos(theta)));
G2 = cos(theta)./(r*(R2+r*cos(theta)));

%hold on
figure(1);
plot(theta, G2, theta, G1, [0,2*pi], [0,0]);
title('Coupling Strength');
xlim([0 2*pi]);
%hold off
legend('R=40','R=80');
xlabel('$\theta$','Interpreter','LaTex')