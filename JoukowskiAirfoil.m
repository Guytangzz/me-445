%-------------------------------------------------------------------------
%   Program for Joukowski Mapping // transformation 
%-------------------------------------------------------------------------
%  Author : Modified by vasishta bhargava 
%-------------------------------------------------------------------------

clear all
close all
clc

%-------------------------------------------
% GRID & INPUT VARIABLES
%------------------------------------------
U = 38
v =U/U;                     % input free stream velocity [m/s] normalized 
alpha = 10;                % angle of attack, deg
alpha = alpha*pi/180;       % angle of attack in radians
sx = 0.0;                 % displacement of circle center in  real axis. // velocity potential
sy = -0.0    ;                 % displacement of circle center in  imaginary axis.  // stream function
s = sx + i*sy;              % resultant displacement in the z plane 
r = 0.35;                   % radius of the cylinder in complex z plane (x+iy)    
 
rho = 1.225;                % fluid density in kg/m3

lambda = r-s;               % offset distance in the complex plane (x+iy)

beta = (alpha);               
k = -2*r*v*sin(beta);        % vortex strength
Gamma = k/(2*pi);           % circulation 

tol = +5e-2;                % geometric grid tolerance for flow visualization 


[x y]= meshgrid(-2:.01:2);   % mesh or grid generation in the circle or complex z plane 
z = x + i*y;


%-------------
% MAIN 
%-------------

w = v * exp(1i*alpha);       % free stream rotation to angle of attack, deg

% grid tolerance check for flow visualization 
for p = 1:length(x)
    for q = 1:length(y)
        if abs(z(p,q)-s) <=  r - tol
            z(p,q) = NaN;
        end
    end
end


%Total complex aerodynamic potential function and grid in airfoil plane 
for p = 1:length(x)
    for q = 1:length(y)
            f(p,q) = w*(z(p,q)) + (v*exp(-1i*alpha)*r^2)./(z(p,q)-s) + 1i*k*log(z(p,q)); 
            J(p,q) = z(p,q)+lambda^2./z(p,q);  % grid in the airfoil plane
    end 
end 


phi = 0:.05:2*pi;   %% JOUKOWSKI MAPPING FUNCTION 
for p = 1:length(phi)
    
        zcirc(p) = r*(cos(phi(p))+1i*sin(phi(p))) + s;
        zair(p) = -zcirc(p)+lambda^2./(-zcirc(p));
end 

% Plotting the solution 
%---------------------------------- 
figure(1)
hold on
v2 = -50:0.1:150;
contourf(real(z),imag(z),imag(f),v2);  
fill(real(zcirc),imag(zcirc),'k')
axis equal
axis(1.5*[-1 1 -1 1])
title('Stream lines around the cylinder.');


figure(2)
hold on
contourf(real(J),imag(J),imag(f),v2); 
fill(real(zair),imag(zair),'k')
axis equal
axis(1.5*[-1 1 -1 1])
title('Stream line contour on mapped airfoil. ');

