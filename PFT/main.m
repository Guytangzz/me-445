clc
clear all

%% Initialization
U_inf = 10; % [m/s] Free-stream velocity
rho = 1.225; % [kg/m^3] Density of air
c = 1; % [m] Chord Length
a = c/4;
alpha = -4;
theta = linspace(0, 2*pi, 100); % Variable to describe the position on the airfoil

%% Classic Joukowski

% Parameters for camber and thickness
hc = 0.04; 
tc = 0.12;

b = tc/3/sqrt(3)*c; 
lambda = hc/2*c;

zeta_0 = -b + 1i*lambda; % Center of circle in the zeta-plane

% Mapping
R = sqrt((a+b)^2 + lambda^2);
zeta_circ = zeta_0 + R.*exp(1i*theta);
z_wing_joukowski = zeta_circ + a^2./zeta_circ;

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

%% Calling complex potential
Gamma = (4*pi*U_inf*R).*sin(alpha + asin(lambda/R));

[x, y] = meshgrid(linspace(-5, 5, 400), linspace(-5, 5, 400));
z = x + 1i*y;

[u, v] = complexPotential(U_inf, alpha, R, Gamma, z, a, zeta_0);

[u_w, v_w] = complexPotential(U_inf, alpha, R, Gamma, z_wing_joukowski, a, zeta_0);

Cp = 1 - sqrt(u_w.^2 + v_w.^2)/U_inf.^2;


%%
figure(1)
contour(x, y, sqrt(u.^2 + v.^2));
xlim([-5, 5]);
ylim([-5, 5]);
grid on;
hold on;
%fill(real(z_wing_joukowski),imag(z_wing_joukowski),'y')
%%
% Create a quiver plot
% figure(1);
% plot([x_up_approx, flip(x_down_approx)], Cp_w, 'b');
% title('Complex Potential Function');
% xlabel('Real Part');
% ylabel('Imaginary Part');
% % xlim([-2, 2]);
% % ylim([-2, 2]);
% % grid on;
% % hold on;
% % fill(real(z_wing_joukowski),imag(z_wing_joukowski),'y')

figure(2);
contour(x, y, sqrt(u.^2 + v.^2), 'b');
title('Complex Potential Function');
xlabel('Real Part');
ylabel('Imaginary Part');
xlim([-5, 5]);
ylim([-5, 5]);
grid on;
hold on;
fill(real(z_wing_joukowski),imag(z_wing_joukowski),'y')

figure(3);
contour(x, y, sqrt(u_s.^2 + v_s.^2), 'b');
title('Complex Potential Function');
xlabel('Real Part');
ylabel('Imaginary Part');
xlim([-5, 5]);
ylim([-5, 5]);
grid on;
hold on;
fill(real(z_wing_joukowski),imag(z_wing_joukowski),'y')

figure(4)
contour(x, y, sqrt((u+u_s).^2 + (v+v_s).^2), 'b');
title('Complex Potential Function');
xlabel('Real Part');
ylabel('Imaginary Part');
xlim([-5, 5]);
ylim([-5, 5]);
grid on;
hold on;
fill(real(z_wing_joukowski),imag(z_wing_joukowski),'y')

figure(5);
contour(x, y, sqrt(u_total.^2 + v_total.^2), 'b');
title('Complex Potential Function with Ground Effect');
xlabel('Real Part');
ylabel('Imaginary Part');
xlim([-5, 5]);
ylim([-5, 5]);
grid on;
hold on;
fill(real(z_wing_joukowski), imag(z_wing_joukowski), 'y');
fill(real(z_ground), imag(z_ground), 'g', 'FaceAlpha', 0.3); % Ground