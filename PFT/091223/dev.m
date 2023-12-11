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

[x, y] = meshgrid(linspace(-2, 2, 400), linspace(-2, 2, 400));
zeta_plane = x + 1i*y;

for ind_a = 1:length(x)
    for ind_b = 1:length(y)
        if abs(zeta_plane(ind_a,ind_b)-zeta_0) <=  R_joukowski
            zeta_plane(ind_a,ind_b) = NaN;
        end
    end
end

z_plane_joukowski = zeta_plane + a^2./zeta_plane;

[x_mir, y_mir] = meshgrid(linspace(-2, 2, 400), linspace(-10, 10, 2000));
zeta_plane_mir = x_mir + 1i*y_mir;

for ind_a = 1:length(zeta_plane_mir(:,1))
    for ind_b = 1:length(zeta_plane_mir(1,:))
        if abs(zeta_plane_mir(ind_a,ind_b)-conj(zeta_0)) <=  R_joukowski
            zeta_plane_mir(ind_a,ind_b) = NaN;
        end
    end
end

z_plane_mir = zeta_plane_mir + a^2./zeta_plane_mir;

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
W_pot = zeros(length(theta), length(alpha_jouk));
u_joukowski = zeros(length(theta), length(alpha_jouk));
v_joukowski = zeros(length(theta), length(alpha_jouk));
V_joukowski = zeros(length(theta), length(alpha_jouk));

Cp_joukowski = zeros(length(theta), length(alpha_jouk));
Cp_up_joukowski = zeros(length(theta)/2, length(alpha_jouk));
Cp_down_joukowski = zeros(length(theta)/2, length(alpha_jouk));

% Coefficients computation
for i = 1:length(alpha_jouk)
    W_pot(:, i) = U_inf*(zeta_circ - zeta_0)*exp(-1i*alpha_jouk(i))+U_inf*R_joukowski^2./(zeta_circ - zeta_0)*exp(1i*alpha_jouk(i))+1i*Gamma_joukowski(i)/(2*pi)*log(zeta_circ-zeta_0);
    W_joukowski(:, i) = (U_inf.*exp(-1i.*alpha_jouk(i)) + 1i.*Gamma_joukowski(i)./(2*pi.*(zeta_circ-zeta_0)) - U_inf.*R_joukowski^2.*exp(1i.*alpha_jouk(i))./((zeta_circ-zeta_0).^2))./(1-a^2./((zeta_circ).^2));

    W_derivative_tot = U_inf*(zeta_circ - zeta_0)*exp(-1i*alpha_jouk(i))+U_inf*R_joukowski^2./(zeta_circ - zeta_0)*exp(1i*alpha_jouk(i))+1i*Gamma_joukowski(i)/(2*pi)*log(zeta_circ-zeta_0) + ...
        U_inf*(zeta_circ - conj(zeta_0))*exp(1i*alpha_jouk(i))+U_inf*R_joukowski^2./(zeta_circ - conj(zeta_0))*exp(-1i*alpha_jouk(i))-1i*Gamma_joukowski(i)/(2*pi)*log(zeta_circ-conj(zeta_0));

    u_joukowski(:, i) = real(W_joukowski(:, i));
    v_joukowski(:, i) = -imag(W_joukowski(:, i));
    V_joukowski(:, i) = sqrt(u_joukowski(:, i).^2 + v_joukowski(:, i).^2);
    Cp_joukowski(:, i) = 1 - V_joukowski(:, i).^2/(U_inf)^2;

    Cp_up_joukowski(:, i) = flip(Cp_joukowski(1:length(theta)/2, i)); 
    Cp_down_joukowski(:, i) = Cp_joukowski(length(theta)/2+1:end, i);

    W_pot_plane(:, :, i) = U_inf*(zeta_plane - zeta_0)*exp(-1i*alpha_jouk(i))+U_inf*R_joukowski^2./(zeta_plane - zeta_0)*exp(1i*alpha_jouk(i))+1i*Gamma_joukowski(i)/(2*pi)*log(zeta_plane-zeta_0);
    u_pot_plane(:, :, i) = real(W_pot_plane(:, :, i));
    v_pot_plane(:, :, i) = imag(W_pot_plane(:, :, i));

    W_pot_plane_mir(:, :, i) = U_inf*(zeta_plane_mir - conj(zeta_0))*exp(1i*alpha_jouk(i))+U_inf*R_joukowski^2./(zeta_plane_mir - conj(zeta_0))*exp(-1i*alpha_jouk(i))-1i*Gamma_joukowski(i)/(2*pi)*log(zeta_plane_mir-conj(zeta_0));
    u_pot_plane_mir(:, :, i) = real(W_pot_plane_mir(:, :, i));
    v_pot_plane_mir(:, :, i) = imag(W_pot_plane_mir(:, :, i));

    W_joukowski_plane(:, :, i) = (U_inf.*exp(-1i.*alpha_jouk(i)) + 1i.*Gamma_joukowski(i)./(2*pi.*(zeta_plane-zeta_0)) - U_inf.*R_joukowski^2.*exp(1i.*alpha_jouk(i))./((zeta_plane-zeta_0).^2))./(1-a^2./(zeta_plane.^2));
    u_joukowski_plane(:, :, i) = real(W_joukowski_plane(:, :, i));
    v_joukowski_plane(:, :, i) = -imag(W_joukowski_plane(:, :, i));
    V_joukowski_plane(:, :, i) = sqrt(u_joukowski_plane(:, :, i).^2 + v_joukowski_plane(:, :, i).^2);

    W_pot_tot(:, :, i) = W_pot_plane(:, :, i) + W_pot_plane_mir(700:1099, :, i);
    u_pot_tot(:, :, i) = real(W_pot_tot(:, :, i));
    v_pot_tot(:, :, i) = imag(W_pot_tot(:, :, i));
end

%% RESULTS OF PFT

%% Cp plot for alpha = -4° joukowski
index_a5_jouk = find(rad2deg(alpha_jouk) == 0);
x_a5 = [x_up_approx, flip(x_down_approx)];

Cp_a5_jouk = [Cp_up_joukowski(:, index_a5_jouk)', flip(Cp_down_joukowski(:, index_a5_jouk)')];

% Create the first square figure (Figure 2)
width_square = 400;  % specify the width for the square figure in pixels
height_square = 400; % specify the height for the square figure in pixels

figure('Position', [100, 100, width_square, height_square]);

contour(real(z_plane_joukowski), imag(z_plane_joukowski), v_pot_plane(:, :, index_a5_jouk), 100)
grid on;
hold on;
fill(real(z_wing_joukowski), imag(z_wing_joukowski), 'y')

title('Joukowski Transformation');
xlabel('Real Axis');
ylabel('Imaginary Axis');

% Create the second rectangular figure with an aspect ratio of 4 (Figure 3)
width_rectangular = 400;   % specify the width for the rectangular figure in pixels
height_rectangular = 1600;  % specify the height for the rectangular figure in pixels

figure('Position', [100, 100, width_rectangular, height_rectangular]);

contour(real(z_plane_mir), imag(z_plane_mir), v_pot_plane_mir(:, :, index_a5_jouk), 400)
grid on;
hold on;
% ylim([-2, 2]);
% xlim([-2, 2]);
fill(real(z_wing_joukowski), -imag(z_wing_joukowski), 'y')

title('Mirror Image Transformation');
xlabel('Real Axis');
ylabel('Imaginary Axis');

% Create the first square figure (Figure 4)
width_square = 400;  % specify the width for the square figure in pixels
height_square = 400; % specify the height for the square figure in pixels

figure('Position', [100, 100, width_square, height_square]);

contour(real(z_plane_joukowski), imag(z_plane_joukowski), v_pot_tot(:, :, index_a5_jouk), 100)
grid on;
hold on;
fill(real(z_wing_joukowski), imag(z_wing_joukowski), 'y')

title('Joukowski Transformation with Mirror Image');
xlabel('Real Axis');
ylabel('Imaginary Axis');
