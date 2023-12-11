% Parameters for NACA 4412
hc = 0.04;  % Maximum camber
pc = 0.4;   % Position of maximum camber
tc = 0.12;  % Thickness-to-chord ratio

% Parameters for Joukowski transformation
a = c/4;  % Half-chord of Joukowski airfoil
xc = -tc/(3*sqrt(3));
yc = hc/2;
R = sqrt((a + tc/2)^2 + yc^2);

% Generate theta values
theta = linspace(0, 2*pi, 100);

% Joukowski transformation
zeta_circ = xc + 1i*yc + R*exp(1i*theta);
zeta_joukowski = zeta_circ + a^2./zeta_circ;

% Plot the Joukowski airfoil in the zeta-plane
figure;
plot(real(zeta_joukowski), imag(zeta_joukowski), '-o');
axis equal;
xlabel('Real part');
ylabel('Imaginary part');
title('Joukowski Airfoil in \zeta-Plane');
grid on;

display(xc)
display(yc)
display(R)