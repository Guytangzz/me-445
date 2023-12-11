% Define parameters
c = 1.0;   % Chord length
V_inf = 1.0;  % Freestream velocity

% Generate airfoil points using Joukowski transform
theta = linspace(0, 2*pi, 100);
z = c * (cos(theta) + 1i * sin(theta));
w = z + c^2 ./ z;

% Calculate velocity components
dw_dz = 1 - c^2 ./ z.^2;
u = real(dw_dz);
v = imag(dw_dz);

% Calculate pressure coefficient
Cp = 1 - (u.^2 + v.^2) / V_inf^2;

% Plot Cp distribution
figure;
plot(real(z/c), Cp);
xlabel('Chord-wise Coordinate (x/c)');
ylabel('Pressure Coefficient (Cp)');
title('Joukowski Airfoil Cp Distribution');
grid on;
