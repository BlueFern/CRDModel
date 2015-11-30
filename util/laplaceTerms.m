% Plot terms of laplace operator and compare magnitude

clear;

theta = 0:0.01:2*pi;
R = 40/(2*pi);
r = 20/(2*pi);

% Term in front of d/dtheta
y1 = -sin(theta)./(r*(R+r*cos(theta)));

% Term in front of d2/dtheta2
y2 = 1/(r^2)*(theta./theta);

% Term in front of d2/dphi2
y3 = 1./((R+r*cos(theta)).^2);

figure(2);
plot(theta, y1, theta, y2, theta, y3);
legend('\partial/\partial\theta', '\partial^2/\partial\theta^2', '\partial^2/\partial\phi^2', 'Interpreter','LaTex');
xlabel('$\theta$','Interpreter','LaTex')
xlim([0 2*pi]);
title('Terms in Laplace operator (R=40)');