clear all
close all
clc

%% Init
v_inf = 10;
v = v_inf/v_inf;
theta = 10;
theta = theta*pi/180;
rho = 1.225;

s_x = -0.0231;
s_y = 0.02;
s = s_x + i*s_y;
r = 0.3106;
lambda = r-s;
beta = (theta);
k = 2*r*v*sin(beta);
Gamma = k/(2*pi);

%% Parameters
w = v * exp(i*theta);

% AERODYNAMIC POTENTIAL
f = @(x,y) w*(x+1i*y) + (v*exp(-1i*theta)*r^2)./(x+1i*y-s) + 1i*k*log(x+1i*y);
% JOUKOWSKI TRANSFORMATION, 
J = @(x, y) x+1i*y+lambda^2./(x+1i*y);
%GRAPHIC - Circle and Joukowski Airfoil
angle = 0:.1:2*pi;
z_circle = r*(cos(angle)+1i*sin(angle)) + s;
z_airfoil = z_circle+lambda^2./z_circle;
% KUTTA JOUKOWSKI THEOREM
L = v_inf*rho*Gamma;
L_str = num2str(L);

contour(real(J(linspace(-5, 5), linspace(-5, 5))), imag(J(linspace(-5, 5), linspace(-5, 5))), imag(f(linspace(-5, 5), linspace(-5, 5))),[-5:.2:5])

