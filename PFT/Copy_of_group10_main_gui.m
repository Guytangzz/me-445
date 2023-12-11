
%%
clc;
clear;
close all;

%%
%%%%%%%%%%%%%%%%%%% POTENTIAL FLOW THEORY %%%%%%%%%%%%%%%%%%%

%% Parameters

U_inf = 10; %[m/s] % Free-stream velocity
rho = 1.225; %[kg/m^3] % Density of air
c=1; %[m] % Chord Length
a = c/4; 
theta=linspace(0,2*pi,100); % Variable to describe the position on the airfoil

%%
%%%%%%%%%%%%%%%%%%%  Classic Joukowski  %%%%%%%%%%%%%%%%%%%

%parameters required to get the good values of max camber and thickness 
hc = 0.04; 
tc = 0.12;

b = tc/3/sqrt(3)*c; 
lambda = hc/2*c; % From question f) of exercice 5) of tutorial on Potential Flow Theory

zeta_0 = -b + 1i*lambda; % Center of circle in the zeta-plane

%Mapping
R_joukowski = sqrt((a+b)^2+lambda^2); % Radius of the circle in the zeta-plane
zeta_circ = zeta_0 + R_joukowski.*exp(1i*theta); % parametrisation of the circle in the zeta-plane
z_wing_joukowski= zeta_circ + a^2./zeta_circ; % Joukowski transformation

%Approximation
m = sqrt(lambda^2+b^2);
eps = m/a;
delta = pi-atan(lambda/b);

x_wing_approx=2*a*cos(theta); % x-values of the Joukowski airfoil in the z-plane
y_wing_approx=2*a*eps.*(cos(delta-theta)-cos(delta)).*sin(theta); % y-values of the Joukowski airfoil in the z-plane

%Obtain x coordinatates form 0 to c & mean camber line 
x_up_approx = flip(x_wing_approx(1:length(theta)/2)+0.5);
x_down_approx = flip(x_wing_approx(end:-1:length(theta)/2+1)+0.5);

y_up_approx = y_wing_approx(1:length(theta)/2);
y_down_approx = y_wing_approx(end:-1:(length(theta)/2+1));

yc_approx = 0.5*(y_up_approx+y_down_approx); %mean camber line 

%% Plotting shapes (AOA 0°)

% Classical Joukowski airfoil

fig_jouk_process = figure; %classic joukowski process and approx
subplot(1,2,1), plot(zeta_circ), axis equal, grid on
xlabel('\xi'), ylabel('\eta'), title('\zeta_{plane}')
subplot(1,2,2), hold all, plot(z_wing_joukowski), plot(x_wing_approx, y_wing_approx)
axis equal, grid on
xlabel('x'), ylabel('y'), title('z_{plane}')
legend('direct Joukowski transform','approximation')
saveas(fig_jouk_process,'fig_jouk_process','png')

fig_jouk_approx_camber = figure; %approx at origin with mean camber line
plot(flip(x_up_approx),y_up_approx, 'blue')
hold on 
plot(flip(x_down_approx),y_down_approx, 'blue')
plot(x_down_approx,yc_approx)
hold off
title('"Equivalent" joukowski airfoil (same thickness, same max camber)')
grid on
axis image

%% Coefficients computation (Joukowski)

alpha_jouk = deg2rad(-6:20); %[rad] % Range of Angle of Attack

Gamma_joukowski = (4*pi*U_inf*R_joukowski).*sin(alpha_jouk + asin(lambda/R_joukowski)); % Kutta Condition;

%pre allocating space for increased running speed
W_joukowski = zeros(length(theta),length(alpha_jouk));
u_joukowski = zeros(length(theta),length(alpha_jouk));
v_joukowski = zeros(length(theta),length(alpha_jouk));
V_joukowski = zeros(length(theta),length(alpha_jouk));

Cp_joukowski = zeros(length(theta),length(alpha_jouk));
Cp_up_joukowski = zeros(length(theta)/2,length(alpha_jouk));
Cp_down_joukowski = zeros(length(theta)/2,length(alpha_jouk));

c_n_ba = zeros(1,length(alpha_jouk)); %ba standing for Basic Aerodynamics
cl_ba = zeros(1,length(alpha_jouk)); 
cl_kj = zeros(1,length(alpha_jouk)); %kj standing for Kutta Joukowski 
cm_le_ba = zeros(1,length(alpha_jouk));
cm_14_joukowski = zeros(1,length(alpha_jouk));

%Coefficients computation
for i = 1:length(alpha_jouk) % For every alpha
    
    W_joukowski(:,i) = (U_inf.*exp(-1i.*alpha_jouk(i))+1i.*Gamma_joukowski(i)./(2*pi.*(zeta_circ-zeta_0))-U_inf.*R_joukowski^2.*exp(1i.*alpha_jouk(i))./((zeta_circ-zeta_0).^2))./(1-a^2./(zeta_circ.^2)); % Dw/Dz (=W)
    u_joukowski(:,i) = real(W_joukowski(:,i)); % x-component of the velocity
    v_joukowski(:,i) = -imag(W_joukowski(:,i));% y-component of the velocity
    V_joukowski(:,i) = sqrt(u_joukowski(:,i).^2+v_joukowski(:,i).^2); % Velocity norm
    Cp_joukowski(:,i) = 1-V_joukowski(:,i).^2/(U_inf)^2; % Pressure coefficients

    Cp_up_joukowski(:,i) = flip(Cp_joukowski(1:length(theta)/2,i)); 
    Cp_down_joukowski(:,i) = Cp_joukowski(length(theta)/2+1:end,i);

    c_n_ba(i) = trapz(x_down_approx,Cp_down_joukowski(:,i))-trapz(x_up_approx,Cp_up_joukowski(:,i)); % normal force coefficients using basic aerodynamics
    cl_ba(i) = cos(alpha_jouk(i))*c_n_ba(i); % for small angles of attack

    L_kj = rho*U_inf*Gamma_joukowski(i); % Lift using the Kutta-Joukowski Theorem
    cl_kj(i) = 2*L_kj/(rho*U_inf^2*c); % with Kutta-Joukouwski lift
    
    
    cm_le_ba(i) = trapz(x_up_approx,Cp_up_joukowski(:,i).*x_up_approx')-trapz(x_down_approx,Cp_down_joukowski(:,i).*x_down_approx'); % Leading-edge moment coefficient using basic aerodynamics equations
    cm_14_joukowski(i) = cm_le_ba(i) + 0.25*c_n_ba(i); % quarter-chord moment coefficient

end

%% RESULTS OF PFT


%% Cp plot for alpha = -4° joukowski

index_a5_jouk = find(rad2deg(alpha_jouk) == -4);
x_a5 = [x_up_approx, flip(x_down_approx)];

Cp_a5_jouk = [Cp_up_joukowski(:,index_a5_jouk)', flip(Cp_down_joukowski(:,index_a5_jouk)')];

fig_cp_jouk = figure;
plot(Cp_a5_jouk,'-ob')
xlim([0,1])
xlabel('x/c')
ylabel('C_p')
legend('Cp_{theo}','location','southeast')
title('Joukowski Mapping')
grid on
set(gca,'TickLabelInterpreter','latex')
set(gcf,'color','w')

saveas(fig_cp_jouk,'fig_cp_jouk','png')

%%
Cp_5 = Cp_joukowski(:, index_a5_jouk);
upper = Cp_5(1:51);
lower = Cp_5(51:end);
plot(x_a5(51:end), -flip(lower))
hold on;
plot(x_a5(1:51), -flip(upper))
legend('lower', 'upper')