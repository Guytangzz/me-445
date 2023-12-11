clc;
clear;
close all;

%% POTENTIAL FLOW THEORY

% Parameters
U_inf = 10; % [m/s] Free-stream velocity
rho = 1.225; % [kg/m^3] Density of air
c = 1; % [m] Chord Length
a = c/4; 
theta = linspace(0, 2*pi, 100); % Variable to describe the position on the airfoil

%% Classic Joukowski

% Parameters for camber and thickness
hc = 0.04; 
tc = 0.12;

b = tc/3/sqrt(3)*c; 
lambda = hc/2*c;

zeta_0 = -b + 1i*lambda; % Center of circle in the zeta-plane

% Mapping
R_joukowski = sqrt((a+b)^2 + lambda^2);
zeta_circ = zeta_0 + R_joukowski.*exp(1i*theta);
z_wing_joukowski = zeta_circ + a^2./zeta_circ;

[x, y] = meshgrid(linspace(-5, 5, 400), linspace(-5, 5, 400));
zeta_plane = x + 1i*y;

% Approximation
m = sqrt(lambda^2 + b^2);
eps = m/a;
delta = pi - atan(lambda/b);

x_wing_approx = 2*a*cos(theta);
y_wing_approx = 2*a*eps.*(cos(delta-theta)-cos(delta)).*sin(theta);

% Obtain x coordinates from 0 to c and mean camber line 
x_up_approx = flip(x_wing_approx(1:length(theta)/2) + 0.5);
x_down_approx = flip(x_wing_approx(end:-1:length(theta)/2+1) + 0.5);

y_up_approx = y_wing_approx(1:length(theta)/2);
y_down_approx = y_wing_approx(end:-1:(length(theta)/2+1));

yc_approx = 0.5*(y_up_approx + y_down_approx); % Mean camber line 

%% Plotting shapes (AOA 0°)

% Classical Joukowski airfoil
fig_jouk_process = figure;
subplot(1,2,1), plot(zeta_circ), axis equal, grid on
xlabel('\xi'), ylabel('\eta'), title('\zeta_{plane}')
subplot(1,2,2), hold all, plot(z_wing_joukowski), plot(x_wing_approx, y_wing_approx)
axis equal, grid on
xlabel('x'), ylabel('y'), title('z_{plane}')
legend('Direct Joukowski transform', 'Approximation')

fig_jouk_approx_camber = figure;
plot(flip(x_up_approx), y_up_approx, 'blue')
hold on 
plot(flip(x_down_approx), y_down_approx, 'blue')
plot(x_down_approx, yc_approx)
hold off
title('"Equivalent" Joukowski airfoil (same thickness, same max camber)')
grid on
axis image

%% Coefficients computation (Joukowski)

alpha_jouk = deg2rad(-6:20); % Range of Angle of Attack

Gamma_joukowski = (4*pi*U_inf*R_joukowski).*sin(alpha_jouk + asin(lambda/R_joukowski));

% Preallocate space
W_joukowski = zeros(length(theta), length(alpha_jouk));
u_joukowski = zeros(length(theta), length(alpha_jouk));
v_joukowski = zeros(length(theta), length(alpha_jouk));
V_joukowski = zeros(length(theta), length(alpha_jouk));

Cp_joukowski = zeros(length(theta), length(alpha_jouk));
Cp_up_joukowski = zeros(length(theta)/2, length(alpha_jouk));
Cp_down_joukowski = zeros(length(theta)/2, length(alpha_jouk));
u_z_plane = zeros(size(zeta_plane));
v_z_plane = zeros(size(zeta_plane));
    
% Coefficients computation
for i = 1:length(alpha_jouk)
    W_joukowski(:, i) = (U_inf.*exp(-1i.*alpha_jouk(i)) + 1i.*Gamma_joukowski(i)./(2*pi.*(zeta_circ-zeta_0)) - U_inf.*R_joukowski^2.*exp(1i.*alpha_jouk(i))./((zeta_circ-zeta_0).^2))./(1-a^2./(zeta_circ.^2));
    u_joukowski(:, i) = real(W_joukowski(:, i));
    v_joukowski(:, i) = -imag(W_joukowski(:, i));
    V_joukowski(:, i) = sqrt(u_joukowski(:, i).^2 + v_joukowski(:, i).^2);
    Cp_joukowski(:, i) = 1 - V_joukowski(:, i).^2/(U_inf)^2;

    Cp_up_joukowski(:, i) = flip(Cp_joukowski(1:length(theta)/2, i)); 
    Cp_down_joukowski(:, i) = Cp_joukowski(length(theta)/2+1:end, i);

    % Use similar coefficients computation as in the original code
    W_z_plane = (U_inf.*exp(-1i.*alpha_jouk(i)) + 1i.*Gamma_joukowski(i)./(2*pi.*(zeta_plane-zeta_0)) - U_inf.*R_joukowski^2.*exp(1i.*alpha_jouk(i))./((zeta_plane-zeta_0).^2))./(1-a^2./(zeta_plane.^2));
    u_z_plane(:, i) = real(W_z_plane);
    v_z_plane(:, i) = -imag(W_z_plane);
end

%% RESULTS OF PFT

%% Cp plot for alpha = -4° joukowski
index_a5_jouk = find(rad2deg(alpha_jouk) == -4);
x_a5 = [x_up_approx, flip(x_down_approx)];

Cp_a5_jouk = [Cp_up_joukowski(:, index_a5_jouk)', flip(Cp_down_joukowski(:, index_a5_jouk)')];

fig_cp_jouk = figure;
plot(x_a5, Cp_a5_jouk, '-ob')
xlim([0, 1])
xlabel('x/c')
ylabel('C_p')
legend('Cp_{theo}', 'Location', 'southeast')
title('Joukowski Mapping')
grid on
set(gca, 'TickLabelInterpreter', 'latex', 'YDir', 'reverse')
set(gcf, 'color', 'w')

figure;
contour(x, y, u_z_plane, 'LineWidth', 1.5);
hold on;
contour(x, y, v_z_plane, 'LineWidth', 1.5);
xlabel('x');
ylabel('y');
title('Velocity Components u and v for z_{plane}');
legend('u_{z_{plane}}', 'v_{z_{plane}}');
axis equal;
grid on;
