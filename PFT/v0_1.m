clc;
clear;
close all;

%% POTENTIAL FLOW THEORY

% Parameters
U_inf = 10; % [m/s] Free-stream velocity
rho = 1.225; % [kg/m^3] Density of air
c = 1; % [m] Chord Length
a = c/4;
h = 0.1; % Height above the ground

theta = linspace(0, 2*pi, 100); % Variable to describe the position on the airfoil

%% Classic Joukowski

% Parameters for camber and thickness
hc = 0.04;
tc = 0.12;

b = tc/3/sqrt(3)*c;
lambda = hc/2*c;

zeta_0 = -b + 1i*lambda; % Center of circle in the zeta-plane

% Mapping with ground effect (virtual image)
R_joukowski = sqrt((a+b)^2 + lambda^2);
zeta_circ = zeta_0 + R_joukowski.*exp(1i*theta);
z_wing_joukowski = zeta_circ + a^2./zeta_circ;

% Image of the airfoil below the ground
zeta_circ_image = zeta_0 - R_joukowski.*exp(1i*theta);
z_wing_joukowski_image = zeta_circ_image + a^2./zeta_circ_image;

% Combine airfoil and image
z_wing_joukowski_combined = [z_wing_joukowski, z_wing_joukowski_image];
theta_combined = [theta, theta];

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
subplot(1,2,2), hold all, plot(z_wing_joukowski_combined), plot(x_wing_approx, y_wing_approx)
axis equal, grid on
xlabel('x'), ylabel('y'), title('z_{plane}')
legend('Direct Joukowski transform', 'Approximation', 'Location', 'southeast')

fig_jouk_approx_camber = figure;
plot(flip(x_up_approx), y_up_approx, 'blue')
hold on 
plot(flip(x_down_approx), y_down_approx, 'blue')
plot(x_down_approx, yc_approx)
hold off
title('"Equivalent" Joukowski airfoil (same thickness, same max camber)')
grid on
axis image

%% Coefficients computation (Joukowski) with ground effect

alpha_jouk = deg2rad(-6:20); % Range of Angle of Attack

Gamma_joukowski = (4*pi*U_inf*R_joukowski).*sin(alpha_jouk + asin(lambda/R_joukowski));

% Preallocate space
W_joukowski_combined = zeros(length(theta_combined), length(alpha_jouk));
u_joukowski_combined = zeros(length(theta_combined), length(alpha_jouk));
v_joukowski_combined = zeros(length(theta_combined), length(alpha_jouk));
V_joukowski_combined = zeros(length(theta_combined), length(alpha_jouk));

Cp_joukowski_combined = zeros(length(theta_combined), length(alpha_jouk));
Cp_up_joukowski_combined = zeros(length(theta_combined)/2, length(alpha_jouk));
Cp_down_joukowski_combined = zeros(length(theta_combined)/2, length(alpha_jouk));

c_n_ba_combined = zeros(1, length(alpha_jouk));
cl_ba_combined = zeros(1, length(alpha_jouk)); 
cl_kj_combined = zeros(1, length(alpha_jouk));
cm_le_ba_combined = zeros(1, length(alpha_jouk));
cm_14_joukowski_combined = zeros(1, length(alpha_jouk));

% Coefficients computation with ground effect
% Coefficients computation with ground effect
for i = 1:length(alpha_jouk)
    % Modified W for ground effect
    W_joukowski_combined(:, i) = (U_inf.*exp(-1i.*alpha_jouk(i)) + ...
        1i.*Gamma_joukowski(i)./(2*pi.*(z_wing_joukowski_combined-zeta_0)) - ...
        U_inf.*R_joukowski^2.*exp(1i.*alpha_jouk(i))./((z_wing_joukowski_combined-zeta_0).^2))./(1-a^2./(z_wing_joukowski_combined.^2));
    
    % Separate into airfoil and image contributions
    W_joukowski(:, i) = W_joukowski_combined(1:length(theta), i);
    W_joukowski_image(:, i) = W_joukowski_combined(length(theta)+1:end, i);
    
    u_joukowski_combined(:, i) = real(W_joukowski_combined(:, i));
    v_joukowski_combined(:, i) = -imag(W_joukowski_combined(:, i));
    V_joukowski_combined(:, i) = sqrt(u_joukowski_combined(:, i).^2 + v_joukowski_combined(:, i).^2);
    Cp_joukowski_combined(:, i) = 1 - V_joukowski_combined(:, i).^2/(U_inf)^2;

    Cp_up_joukowski_combined(:, i) = flip(Cp_joukowski_combined(1:length(theta_combined)/2, i)); 
    Cp_down_joukowski_combined(:, i) = Cp_joukowski_combined(length(theta_combined)/2+1:end, i);

    % Use x_wing_approx instead of x_up_approx and x_down_approx
    c_n_ba_combined(i) = trapz(x_wing_approx, Cp_down_joukowski_combined(:, i)) - trapz(x_wing_approx, Cp_up_joukowski_combined(:, i));
    cl_ba_combined(i) = cos(alpha_jouk(i)) * c_n_ba_combined(i);

    L_kj_combined = rho * U_inf * Gamma_joukowski(i);
    cl_kj_combined(i) = 2 * L_kj_combined / (rho * U_inf^2 * c);
    
    cm_le_ba_combined(i) = trapz(x_wing_approx, Cp_up_joukowski_combined(:, i) .* x_wing_approx') - trapz(x_wing_approx, Cp_down_joukowski_combined(:, i) .* x_wing_approx');
    cm_14_joukowski_combined(i) = cm_le_ba_combined(i) + 0.25 * c_n_ba_combined(i);
end

%% RESULTS OF PFT

%% Cp plot for alpha = -4° joukowski
index_a5_jouk = find(rad2deg(alpha_jouk) == -4);
x_a5_combined = [x_up_approx, flip(x_down_approx)];

Cp_a5_jouk_combined = [Cp_up_joukowski_combined(:, index_a5_jouk)', flip(Cp_down_joukowski_combined(:, index_a5_jouk)')];

fig_cp_jouk_combined = figure;
plot(x_a5_combined, Cp_a5_jouk_combined, '-ob')
xlim([0, 1])
xlabel('x/c')
ylabel('C_p')
legend('Cp_{theo}', 'Location', 'southeast')
title('Joukowski Mapping with Ground Effect')
grid on
set(gca, 'TickLabelInterpreter', 'latex', 'YDir', 'reverse')
set(gcf, 'color', 'w')
