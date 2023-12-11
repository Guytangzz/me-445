
%%
clc;
clear;
close all;

%% Loading data

% Load (x/c; y/c) coordinates of the airfoil geometry from external text
% file
fileName = 'S809.txt';

airfoil_geometry = readmatrix(fileName);

% upper and lower curves of the airfoil
upper_curve=sortrows(airfoil_geometry(1:33,:),1);
lower_curve=sortrows(airfoil_geometry(34:66,:),1);

% Redefine x/c and y/c for upper and lower curves
x_c_TAT=linspace(0,1,101);
y_c_l_TAT=interp1(lower_curve(:,1),lower_curve(:,2),x_c_TAT,'linear','extrap');
y_c_u_TAT=interp1(upper_curve(:,1),upper_curve(:,2),x_c_TAT,'linear','extrap');

% camber line  
ycamber_c_TAT=(y_c_u_TAT+y_c_l_TAT)/2;

% plot airfoil geometry and camberline
fig_geom = figure;
hold on
plot(x_c_TAT,ycamber_c_TAT,'b-','LineWidth',1.5);
plot(lower_curve(:,1),lower_curve(:,2),'k-','LineWidth',1.5)
plot(upper_curve(:,1),upper_curve(:,2),'k-','LineWidth',1.5)
hold off
title('S809 AirfoilTT','Interpreter','latex')
xlabel('$x/c$','Interpreter','latex')
ylabel('$y/c$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'color','w');
axis image

saveas(fig_geom, 'fig_geom','png')

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

%%

%%%%%%%%%%%%%%%%%%%  Alternative Mapping  %%%%%%%%%%%%%%%%%%%
%Alternative mappping of the S809 airfoil (prefect shape)

%data extracted from Xudong
A = [1.47486 -0.55641 -0.07436];
B = [0.04733 -0.65773 -0.05212];

C = [1-cos(theta); (1-cos(theta)).^2; (1-cos(theta)).^3];
S = [sin(theta); sin(theta).^2; sin(theta).^3];

phi = A*C + B*S;

%Mapping
R_map = a; % Radius of the circle in the zeta-plane
zeta2 = R_map.*exp(1i*theta); % zeta"  original circle
zeta1 = zeta2.*exp(phi); %zeta'  near circle "shrinked" with the phi transformation
z_wing = zeta1 + a^2./zeta1; % Joukowski transformation

%Obtaining x coordinatates form 0 to c and the mean camber line
x_wing_map = real(z_wing);
y_wing_map =imag(z_wing);

x_up_map = flip(x_wing_map(1:length(theta)/2)+c/2);
x_down_map = flip(x_wing_map(end:-1:length(theta)/2+1)+c/2);

y_up_map = y_wing_map(1:length(theta)/2);
y_down_map = y_wing_map(end:-1:(length(theta)/2+1));

yc_map = 0.5*(y_up_map+y_down_map); %mean camber line


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

% Alternative mapping

fig_map_process = figure; % alternative mapping process
subplot(1,3,1), plot(zeta2), axis equal, grid on
xlabel('\xi'), ylabel('\eta'), title('\zeta"_{plane}')
subplot(1,3,2), plot(zeta1), axis equal, grid on
xlabel('\xi'), ylabel('\eta'), title('\zeta´_{plane}')
subplot(1,3,3), hold all, plot(z_wing)
axis equal, grid on
xlabel('x'), ylabel('y'), title('z_{plane}')
saveas(fig_map_process,'fig_map_process','png')

fig_map_camber = figure; %alternative mapping at origin with mean camber line
plot(flip(x_up_map),y_up_map, 'blue')
hold on 
plot(flip(x_down_map),y_down_map, 'blue')
plot(x_down_map,yc_map)
hold off
title('S809 with the alternative mapping')
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

%% Coefficients computation ("Perfect maping")

%AOA 
alpha1 = deg2rad(linspace(-6,5,42));
alpha2 = deg2rad(linspace(5.01,20,58));
alpha3= [alpha1,alpha2]; % "new" alpha in order to have 100 elements (as phi) and still have an exact value for alpha=5

%Circulation
Gamma_map = 2*1i*pi*U_inf/a .* exp(-1i.*alpha3 - phi).*(-a^2 + R_map^2.*exp(2.*1i.*alpha3 + 2.*phi)); % Kutta Condition;
Gamma_r = real(Gamma_map); %taking the real part of the circulation 
Gamma_i = imag(Gamma_map); %imaginary part not used

%Lift
L_map = -rho.*Gamma_r.*U_inf; % taking real part of the circulation % Kutta-Joukowski Theorem
%Lift coeff 
Cl_map = L_map./(0.5*rho*c*U_inf^2); 

%Applying AOA to the maping
zeta = zeta1.*exp(1i.*alpha3); 

%Pressure coefficient Cp computation, and other coefficients

%pre allocating space for uncreased running speed
W_map = zeros(length(theta),length(alpha3));
u_map = zeros(length(theta),length(alpha3));
v_map = zeros(length(theta),length(alpha3));
V_map = zeros(length(theta),length(alpha3));
Cp_map = zeros(length(theta),length(alpha3));
Cp_up_map = zeros(length(theta)/2,length(alpha3));
Cp_down_map = zeros(length(theta)/2,length(alpha3));

c_n_map = zeros(1,length(alpha3));
cm_le_map = zeros(1,length(alpha3));
cm_14_map = zeros(1,length(alpha3));


for i = 1:length(alpha3) % for all angles of attack
    W_map(:,i) = (U_inf.*(1-R_map^2./(zeta2.^2))-1i.*Gamma_r(i)./(2.*pi.*zeta2)).*exp(-phi).*exp(-1i.*alpha3(i))./(1-a^2./zeta(i).^2); % Dw/Dz = W
    u_map(:,i) = real(W_map(:,i));
    v_map(:,i) = -imag(W_map(:,i));
    V_map(:,i) = sqrt(u_map(:,i).^2+v_map(:,i).^2); % Velocity norm

    Cp_map(:,i) = 1-V_map(:,i).^2/(U_inf)^2; % Pressure coefficients
    Cp_up_map(:,i) = flip(Cp_map(1:length(theta)/2,i)); %extrados
    Cp_down_map(:,i) = Cp_map(length(theta)/2+1:end,i); %intrados

    c_n_map(i) = trapz(x_down_approx,Cp_down_map(:,i))-trapz(x_up_approx,Cp_up_map(:,i)); % normal force coefficient using basic aerodynamics equations

    cm_le_map(i) = trapz(x_up_approx,Cp_up_map(:,i).*x_up_approx')-trapz(x_down_approx,Cp_down_map(:,i).*x_down_approx'); % Leading-edge moment coefficient using basic aerodynamics equations
    cm_14_map(i) = cm_le_map(i) + 0.25*c_n_map(i); % quarter-chord moment coefficient by shifting of the point of reference
end



%% RESULTS OF PFT


%% Cp plot for alpha = 5° alternative/perfect mapping

load Cp_lower_alpha_5.mat
load Cp_upper_alpha_5.mat % Loading of Experimental values

Cp_exp = [Cp_upper_alpha_5(:,2)',flip(Cp_lower_alpha_5(:,2)')];
x_exp = [Cp_upper_alpha_5(:,1)',flip(Cp_lower_alpha_5(:,1)')];

index_a5 = find(rad2deg(alpha3) == 5);

Cp_a5_map = [Cp_up_map(:,index_a5)', flip(Cp_down_map(:,index_a5)')];
x_a5 = [x_up_approx, flip(x_down_approx)];

fig_cp_map = figure ;
plot(x_a5,Cp_a5_map,'-ob')
hold on
plot(x_exp,Cp_exp,'--r')
xlim([0,1])
xlabel('x/c'), ylabel('C_p'), legend('Cp_{theo}','Cp_{exp}','location','southeast'), title('Alternative Mapping')
grid on
set(gca,'TickLabelInterpreter','latex')
set(gcf,'color','w')
saveas(fig_cp_map,'fig_cp_map','png')

%% Cp plot for alpha = 5° joukowski

index_a5_jouk = find(rad2deg(alpha_jouk) == 5);

Cp_a5_jouk = [Cp_up_joukowski(:,index_a5_jouk)', flip(Cp_down_joukowski(:,index_a5_jouk)')];

fig_cp_jouk = figure;
plot(x_a5,Cp_a5_jouk,'-ob')
hold on
plot(x_exp,Cp_exp,'--r')
xlim([0,1])
xlabel('x/c')
ylabel('C_p')
legend('Cp_{theo}','Cp_{exp}','location','southeast')
title('Joukowski Mapping')
grid on
set(gca,'TickLabelInterpreter','latex')
set(gcf,'color','w')

saveas(fig_cp_jouk,'fig_cp_jouk','png')

%% Aerodynamic center joukowski

d_alpha = rad2deg(alpha_jouk(2)-alpha_jouk(1)); % increment of angle of attack

% pre-allocating space 
m_0 = [];
a_0 = [];

for j = 2:(length(alpha_jouk)-1)
    m_0(j) = (cm_14_joukowski(j+1)-cm_14_joukowski(j-1))/(2*d_alpha); % computing numerically the derivative of the leading-edge moment coefficient with respect to the angle of attack
    a_0(j) = (c_n_ba(j+1)-c_n_ba(j-1))/(2*d_alpha); % computing numerically the derivative of the normal force coefficient with respect to the angle of attack
end
m_0(1) = (cm_14_joukowski(2)-cm_14_joukowski(1))/d_alpha; 
m_0(length(alpha_jouk)) = (cm_14_joukowski(end)-cm_14_joukowski(end-1))/d_alpha;
a_0(1) = (c_n_ba(2)-c_n_ba(1))/d_alpha;
a_0(length(alpha_jouk)) = (c_n_ba(end)-c_n_ba(end-1))/d_alpha;

A_0 = mean(a_0(1:16)); % assuming the normal force coefficient evolves linearly with the angle of attack
M_0 = mean(m_0(1:16)); % assuming the leading-edge moment coefficient evolves linearly with the angle of attack

a_c = -M_0/A_0+1/4; % computation of the aerodynamic center

%% Leading-Edge and Quarter Chord Moments joukowski + Moment at aerodynamic center (should be constant with respect to alpha) 

cm_ac = cm_14_joukowski - (0.25-a_c).*c_n_ba; % computation of the aerodynamic center moment coefficient (should be constant)

fig_cm_all = figure;
plot(rad2deg(alpha_jouk),cm_le_ba,'LineWidth',1.5)
hold on
plot(rad2deg(alpha_jouk),cm_14_joukowski,'LineWidth',1.5)
hold on
plot(rad2deg(alpha_jouk), cm_ac,'LineWidth',1.5)
grid on
legend('C_{m}_{LE}','C_{m}_{1/4}','c_{m_{AC}}')
xlabel('\alpha [deg]')
ylabel('c_m')
xlim([-6,20])
set(gca,'TickLabelInterpreter','latex')
set(gcf,'color','w');


%% Center of pressure joukowski

x_cp = -cm_le_ba./c_n_ba*c; % computation of the location of the center of pressure

fig_xcp_all = figure;
plot(rad2deg(alpha_jouk),x_cp)
grid on
legend('x_{Cp} [m]','location','southeast')

% Dividing by approximately 0 when c_n tends to 0. This explains the sudden
% peak


%% 
%%%%%%%%%%%%%%%% THIN AIRFOIL THEORY %%%%%%%%%%%%%%%%%

%if needed the plots with only the TAT method are commented (%) just below,
%but the combined graphs used in the report are plotted after this section 


%% Fourier coefficients dy_c/dx= A_0 + sum(A_n*cos(n*theta)) for classic S809 airfoil

% Derivative of the camber line along chord length
dy_c_dx_c_TAT=gradient(ycamber_c_TAT,x_c_TAT);

% Change of variable x/c=(1-cos(theta))/2
theta_TAT=acos(1-2*x_c_TAT); % rad

% calculation of first Fourier coefficients
A_0_TAT=1/pi*trapz(theta_TAT,dy_c_dx_c_TAT);
A_1_TAT=2/pi*trapz(theta_TAT,dy_c_dx_c_TAT.*cos(theta_TAT));
A_2_TAT=2/pi*trapz(theta_TAT,dy_c_dx_c_TAT.*cos(2*theta_TAT));

%% Lift coefficient (C_l) vs angle of attack (alpha)
 
% Range of angle of attack
alpha_TAT = -5:0.5:20; % degrees

% Experimental lift coefficient vs angle of attack from paper
load Cl_alpha.mat
Cl_exp_TAT = interp1(Cl_vs_AoA(:,1),Cl_vs_AoA(:,2),alpha_TAT);

% Cl according to Thin Airfoil Theory
Cl_TAT = 2*pi*alpha_TAT*pi/180 + pi*(A_1_TAT-2*A_0_TAT);

% Plot theoretical and experimental lift coefficient vs angle of attack
%figure
%hold on
%plot(alpha_TAT,Cl_TAT,'b-','LineWidth',1.5)
%plot(alpha_TAT,Cl_exp_TAT,'g-','LineWidth',1.5)
% hold off
% legend('$Cl_{theo}$','$Cl_{exp}$','interpreter','latex')
% grid on
% title('Lift coefficient ($C_l$) in function of the angle of attack ($\alpha$) - Classic Airfoil','Interpreter','latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$C_l$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');

% Pitching moment coefficient at leading edge and aerodynamic centre 
% (ac=1/4 for Thin Airfoil Theory) and centre of pressure

Cm_LE_TAT=-pi/2*(alpha_TAT*pi/180-A_0_TAT+A_1_TAT-A_2_TAT/2);
Cm_ac_TAT=-pi/4*(A_1_TAT-A_2_TAT)*ones(length(alpha_TAT)); % constant value

xCP_TAT=1/4+pi/4./Cl_TAT*(A_1_TAT-A_2_TAT);

% figure
% hold on
% plot(alpha_TAT,Cm_LE_TAT,'b-','LineWidth',1.5)
% hold off
% grid on
% title('Pitching moment coefficient at the leading edge ($C_{m,LE}$) in function of the angle of attack ($\alpha$) - Classic Airfoil','Interpreter','latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$C_{m,LE}$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');
% 
% figure
% hold on
% plot(alpha_TAT,Cm_ac_TAT,'b-','LineWidth',1.5)
% hold off
% grid on
% title('Pitching moment coefficient at the aerodynamic centre ($C_{m,ac}$) in function of the angle of attack ($\alpha$) - Classic Airfoil','Interpreter','latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$C_{m,ac}$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');
% 
% figure
% hold on
% plot(alpha_TAT,xCP_TAT,'b-','LineWidth',1.5)
% hold off
% grid on
% title('Position of the centre of pressure ($x_{CP}/c$) in function of the angle of attack ($\alpha$) - Classic airfoil','Interpreter','latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$x_{CP}/c$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');


%% Airfoil with traditional flap

%flap deflection angle
alpha_F_TAT=-5:2.5:10; % degrees

% index of the position of the flap (hinge position at x/c=0.9)
index_flap_TAT=find(x_c_TAT>=0.9,1);

% New Fourier coefficients and lift coefficient for the different flap deflection
for j = 1:length(alpha_F_TAT)
    A_0_F_TAT(j)=A_0_TAT-1/pi*alpha_F_TAT(j)*pi/180*(theta_TAT(end)-theta_TAT(index_flap_TAT));
    A_1_F_TAT(j)=A_1_TAT+2/pi*trapz(theta_TAT(index_flap_TAT:end),-alpha_F_TAT(j)*pi/180.*cos(theta_TAT(index_flap_TAT:end)));
    A_2_F_TAT(j)=A_2_TAT+2/pi*trapz(theta_TAT(index_flap_TAT:end),-alpha_F_TAT(j)*pi/180.*cos(2*theta_TAT(index_flap_TAT:end)));

    % Lift coefficient vs AoA for different flap deflection
    Cl_TAT_F(j,:) = 2*pi*alpha_TAT*pi/180 + pi*(A_1_F_TAT(j)-2*A_0_F_TAT(j));
end

% Load experimental data (C_l vs flap deflection for AoA = 0,
% 2.5 and 5 degrees)
load Cl_trad_AoA_0_vs_alpha_F.mat
Cl_exp_trad_AoA_0_TAT = interp1(Cl_trad_AoA_0_vs_alpha_f(:,1),Cl_trad_AoA_0_vs_alpha_f(:,2),alpha_F_TAT);
load Cl_trad_AoA_2.5_vs_alpha_F.mat
Cl_exp_trad_AoA_2_5_TAT = interp1(Cl_trad_AoA_20x2E5_vs_alpha_F(:,1),Cl_trad_AoA_20x2E5_vs_alpha_F(:,2),alpha_F_TAT);
load Cl_trad_AoA_5_vs_alpha_F.mat
Cl_exp_trad_AoA_5_TAT = interp1(Cl_trad_AoA_5_vs_alpha_F(:,1),Cl_trad_AoA_5_vs_alpha_F(:,2),alpha_F_TAT);

% Plot theoritical and experimental lift coefficient vs flap deflection for
% AoA = 0, 2.5 and 5 degrees
% figure
% hold on
% plot(alpha_F_TAT,Cl_TAT_F(:,find(alpha_TAT==0)),'b-*','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_TAT_F(:,find(alpha_TAT==2.5)),'b-o','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_TAT_F(:,find(alpha_TAT==5)),'b-s','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_exp_trad_AoA_0_TAT,'g-*','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_exp_trad_AoA_2_5_TAT,'g-o','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_exp_trad_AoA_5_TAT,'g-s','LineWidth',1.5);
% hold off
% grid on
% title(' Lift coefficient ($C_l$) in function of the flap deflection ($\alpha_F$) - Traditional flap','Interpreter','latex')
% xlabel('$\alpha_F [^{\circ}]$','Interpreter','latex')
% ylabel('$C_l$','Interpreter','latex')
% legend('Theo $\alpha=0^\circ$','Theo $\alpha=2.5^\circ$','Theo $\alpha=5^\circ$','Exp $\alpha=0^\circ$','Exp $\alpha=2.5^\circ$','Exp $\alpha=5^\circ$','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');

% Pitching moment coefficient at leading edge and aerodynamice centre,
% and centre of pressure for different flap defection

for j=1:length(alpha_F_TAT)
    Cm_LE_TAT_F(j,:)=-pi/2*(alpha_TAT*pi/180-A_0_F_TAT(j)+A_1_F_TAT(j)-A_2_F_TAT(j)/2);
    Cm_ac_TAT_F(j,:)=-pi/4*(A_1_F_TAT(j)-A_2_F_TAT(j))*ones(1,length(alpha_TAT)); 

    xCP_TAT_F(j,:)=1/4+pi/4./Cl_TAT_F(j,:)*(A_1_F_TAT(j)-A_2_F_TAT(j));
end

% figure
% hold on
% for j=1:length(alpha_F_TAT)
% plot(alpha_TAT,Cm_LE_TAT_F(j,:),'LineWidth',1.5)
% end
% hold off
% grid on
% title('Pitching moment coefficient at the leading edge ($C_{m,LE}$) in function of the angle of attack ($\alpha$) - Traditional flap','Interpreter','latex')
% legend('$\alpha_F=-5^\circ$','$\alpha_F=-2.5^\circ$','$\alpha_F=0^\circ$','$\alpha_F=2.5^\circ$','$\alpha_F=5^\circ$','$\alpha_F=7.5^\circ$','$\alpha_F=10^\circ$','Interpreter', 'latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$C_{m,LE}$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');
% 
% figure
% hold on
% for j=1:length(alpha_F_TAT)
% plot(alpha_TAT,Cm_ac_TAT_F(j,:),'LineWidth',1.5)
% end
% hold off
% grid on
% title('Pitching moment coefficient at the aerodynamic centre ($C_{m,ac}$) in function of the angle of attack ($\alpha$) - Traditional flap','Interpreter','latex')
% legend('$\alpha_F=-5^\circ$','$\alpha_F=-2.5^\circ$','$\alpha_F=0^\circ$','$\alpha_F=2.5^\circ$','$\alpha_F=5^\circ$','$\alpha_F=7.5^\circ$','$\alpha_F=10^\circ$','Interpreter', 'latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$C_{m,ac}$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');
% 
% figure
% hold on
% for j=1:length(alpha_F_TAT)
% plot(alpha_TAT,xCP_TAT_F(j,:),'LineWidth',1.5)
% end
% hold off
% grid on
% title('Position of the centre of pressure ($x_{CP}/c$) in function of the angle of attack ($\alpha$)- Traditional flap','Interpreter','latex')
% legend('$\alpha_F=-5^\circ$','$\alpha_F=-2.5^\circ$','$\alpha_F=0^\circ$','$\alpha_F=2.5^\circ$','$\alpha_F=5^\circ$','$\alpha_F=7.5^\circ$','$\alpha_F=10^\circ$','Interpreter', 'latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$x_{CP}/c$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');


%% Airfoil with gurney flap (proposed flap in the paper)


% We consider that the gurney effect is equivalent to an additionnal flap
% deflection angle. The gurney has a height h/c = 0.01 and a distance from
% the hinge of l/c=0.08. We add an additional 0.8 factor since the gurney
% is installed at a 0.8 factor of the flap length

delta_alpha_g_TAT=atand(0.8*0.01/0.08); % degrees

% New Fourier coefficients and lift coefficient for the proposed airfoil  
% with a gurney flap for different flap deflection

for j = 1:length(alpha_F_TAT)
    A_0_Fg_TAT(j)=A_0_TAT-1/pi*(alpha_F_TAT(j)+delta_alpha_g_TAT)*pi/180*(theta_TAT(end)-theta_TAT(index_flap_TAT));
    A_1_Fg_TAT(j)=A_1_TAT+2/pi*trapz(theta_TAT(index_flap_TAT:end),-(alpha_F_TAT(j)+delta_alpha_g_TAT)*pi/180.*cos(theta_TAT(index_flap_TAT:end)));
    A_2_Fg_TAT(j)=A_2_TAT+2/pi*trapz(theta_TAT(index_flap_TAT:end),-(alpha_F_TAT(j)+delta_alpha_g_TAT)*pi/180.*cos(2*theta_TAT(index_flap_TAT:end)));

    % lift coefficient with gurney flap vs AoA for different flap
    % deflection
    Cl_TAT_Fg(j,:) = 2*pi*alpha_TAT*pi/180 + pi*(A_1_Fg_TAT(j)-2*A_0_Fg_TAT(j));
end

% Load experimental data (C_l vs flap deflection for AoA = 0,
% 2.5 and 5 degrees)
load Cl_prop_AoA_0_vs_alpha_f.mat
Cl_exp_prop_AoA_0_TAT= interp1(Cl_prop_AoA_0_vs_alpha_f(:,1),Cl_prop_AoA_0_vs_alpha_f(:,2),alpha_F_TAT);
load Cl_prop_AoA_2_5_vs_alpha_f.mat
Cl_exp_prop_AoA_2_5_TAT = interp1(Cl_prop_AoA_2_5_vs_alpha_f(:,1),Cl_prop_AoA_2_5_vs_alpha_f(:,2),alpha_F_TAT);
load Cl_prop_AoA_5_vs_alpha_f.mat
Cl_exp_prop_AoA_5_TAT = interp1(Cl_prop_AoA_5_vs_alpha_f(:,1),Cl_prop_AoA_5_vs_alpha_f(:,2),alpha_F_TAT);


% figure
% hold on
% plot(alpha_F_TAT,Cl_TAT_Fg(:,find(alpha_TAT==0)),'b-*','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_TAT_Fg(:,find(alpha_TAT==2.5)),'b-o','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_TAT_Fg(:,find(alpha_TAT==5)),'b-s','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_exp_prop_AoA_0_TAT,'g-*','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_exp_prop_AoA_2_5_TAT,'g-o','LineWidth',1.5);
% plot(alpha_F_TAT,Cl_exp_prop_AoA_5_TAT,'g-s','LineWidth',1.5);
% hold off
% grid on
% title(' Lift coefficient ($C_l$) in function of the flap deflection ($\alpha_F$) - Gurney flap','Interpreter','latex')
% legend('Theo $\alpha=0^\circ$','Theo $\alpha=2.5^\circ$','Theo $\alpha=5^\circ$','Exp $\alpha=0^\circ$','Exp $\alpha=2.5^\circ$','Exp $\alpha=5^\circ$','interpreter','latex')
% xlabel('$\alpha_F [^{\circ}]$','Interpreter','latex')
% ylabel('$C_l$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');

% Pitching moment coefficient at leading edge and centre of pressure for
% different flap defection

for j=1:length(alpha_F_TAT)
    Cm_LE_TAT_Fg(j,:)=-pi/2*(alpha_TAT*pi/180-A_0_Fg_TAT(j)+A_1_Fg_TAT(j)-A_2_Fg_TAT(j)/2);
    Cm_ac_TAT_Fg(j,:)=-pi/4*(A_1_Fg_TAT(j)-A_2_Fg_TAT(j))*ones(1,length(alpha_TAT)); 

    xCP_TAT_Fg(j,:)=1/4+pi/4./Cl_TAT_Fg(j,:)*(A_1_Fg_TAT(j)-A_2_Fg_TAT(j));
end

% figure
% hold on
% for j=1:length(alpha_F_TAT)
% plot(alpha_TAT,Cm_LE_TAT_Fg(j,:),'LineWidth',1.5)
% end
% hold off
% grid on
% title('Pitching moment coefficient at the leading edge ($C_{m,LE}$) in function of the angle of attack ($\alpha$) - Gurney flap','Interpreter','latex')
% legend('$\alpha_F=-5^\circ$','$\alpha_F=-2.5^\circ$','$\alpha_F=0^\circ$','$\alpha_F=2.5^\circ$','$\alpha_F=5^\circ$','$\alpha_F=7.5^\circ$','$\alpha_F=10^\circ$','Interpreter', 'latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$C_{m,LE}$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');

% figure
% hold on
% for j=1:length(alpha_F_TAT)
% plot(alpha_TAT,Cm_ac_TAT_Fg(j,:),'LineWidth',1.5)
% end
% hold off
% grid on
% title('Pitching moment coefficient at the aerodynamic centre ($C_{m,ac}$) in function of the angle of attack ($\alpha$) - Gurney flap','Interpreter','latex')
% legend('$\alpha_F=-5^\circ$','$\alpha_F=-2.5^\circ$','$\alpha_F=0^\circ$','$\alpha_F=2.5^\circ$','$\alpha_F=5^\circ$','$\alpha_F=7.5^\circ$','$\alpha_F=10^\circ$','Interpreter', 'latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$C_{m,ac}$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');

% figure
% hold on
% for j=1:length(alpha_F_TAT)
% plot(alpha_TAT,xCP_TAT_Fg(j,:),'LineWidth',1.5)
% end
% hold off
% grid on
% title('Position of the centre of pressure ($x_{CP}/c$) in function of the angle of attack ($\alpha$) - Gurney flap','Interpreter','latex')
% legend('$\alpha_F=-5^\circ$','$\alpha_F=-2.5^\circ$','$\alpha_F=0^\circ$','$\alpha_F=2.5^\circ$','$\alpha_F=5^\circ$','$\alpha_F=7.5^\circ$','$\alpha_F=10^\circ$','Interpreter', 'latex')
% xlabel('$\alpha$ [$^{\circ}$]','Interpreter','latex')
% ylabel('$x_{CP}/c$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gcf,'color','w');




%% PLOTTING FIGURES COMBINING THE 2 METHODS

%if needed the plots with only the TAT method are commented (%) before

%% Cl_alpha plot BA/KJ/alternative mapping/TAT vs Experimental

load Cl_alpha.mat % Loading of experimental result
alpha_exp = Cl_vs_AoA(:,1);
cl_exp = Cl_vs_AoA(:,2);
cl_exp_PFT = interp1(Cl_vs_AoA(:,1),Cl_vs_AoA(:,2),rad2deg(alpha3));

fig_cl_all = figure;
plot(rad2deg(alpha_jouk),cl_ba,'-',rad2deg(alpha_jouk),cl_kj,'-',rad2deg(alpha3),Cl_map,'-','Linewidth',1) %plotting the 3 different curves we computed before
hold on
plot(alpha_TAT,Cl_TAT,'-','Linewidth',1)
hold on
plot(alpha_exp,cl_exp,'r--','Linewidth',1) %exp results
legend('Cl_{BA}','Cl_{KJ}','Cl_{AM}','Cl_{TAT}','Cl_{EXP}','location', 'northwest')
grid on
xlabel('\alpha [deg]')
ylabel('c_l')
set(gcf,'color','w');
saveas(fig_cl_all,'fig_cl_all','png')

%% Cl_alpha plot TAT flap 

fig_cl_flap = figure;
hold on
plot(alpha_F_TAT,Cl_TAT_Fg(:,find(alpha_TAT==0)),'b-*','LineWidth',1);
plot(alpha_F_TAT,Cl_TAT_Fg(:,find(alpha_TAT==2.5)),'b-o','LineWidth',1);
plot(alpha_F_TAT,Cl_TAT_Fg(:,find(alpha_TAT==5)),'b-s','LineWidth',1);
plot(alpha_F_TAT,Cl_exp_prop_AoA_0_TAT,'r-*','LineWidth',1);
plot(alpha_F_TAT,Cl_exp_prop_AoA_2_5_TAT,'r-o','LineWidth',1);
plot(alpha_F_TAT,Cl_exp_prop_AoA_5_TAT,'r-s','LineWidth',1);
plot(alpha_F_TAT,Cl_TAT_F(:,find(alpha_TAT==0)),'k-*','LineWidth',1);
plot(alpha_F_TAT,Cl_TAT_F(:,find(alpha_TAT==2.5)),'k-o','LineWidth',1);
plot(alpha_F_TAT,Cl_TAT_F(:,find(alpha_TAT==5)),'k-s','LineWidth',1);
hold off
grid on
legend('Theo gurney \alpha=0\circ','Theo gurney \alpha=2.5\circ','Theo gurney \alpha=5\circ','Exp gurney \alpha=0\circ','Exp gurney \alpha=2.5\circ','Exp gurney\alpha=5\circ','Theo trad \alpha=0\circ','Theo trad \alpha=2.5\circ','Theo trad \alpha=5\circ')
legend('location','northwest')
xlabel('\alpha_F [\circ]')
ylabel('C_l')
set(gcf,'color','w');
saveas(fig_cl_flap,'fig_cl_flap', 'png')

%% Cm,LE trad flap vs prop flap (gurney)

fig_cmle_flap_gurney = figure;
hold on
plot(alpha_TAT,Cm_LE_TAT_Fg(1,:),'r-','LineWidth',1)
plot(alpha_TAT,Cm_LE_TAT_F(1,:),'b-','LineWidth',1)
plot(alpha_TAT,Cm_LE_TAT_Fg(4,:),'r--','LineWidth',1)
plot(alpha_TAT,Cm_LE_TAT_F(4,:),'b--','LineWidth',1)
plot(alpha_TAT,Cm_LE_TAT_Fg(7,:),'r-.','LineWidth',1)
plot(alpha_TAT,Cm_LE_TAT_F(7,:),'b-.','LineWidth',1)
hold off
grid on
legend('Gurney \alpha_F=-5\circ','Trad \alpha_F=-5\circ','Gurney \alpha_F=2.5\circ','Trad \alpha_F=2.5\circ','Gurney \alpha_F=7.5\circ','Trad \alpha_F=7.5\circ')
xlabel('\alpha [deg]')
ylabel('c_{m,LE}')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'color','w');
saveas(fig_cmle_flap_gurney, 'fig_cmle_flap_gurney', 'png')

%% Cm,ac trad flap vs prop flap (gurney)

fig_cmac_flap_gurney = figure;
hold on
plot(alpha_TAT,Cm_ac_TAT_Fg(1,:),'r-','LineWidth',1)
plot(alpha_TAT,Cm_ac_TAT_F(1,:),'b-','LineWidth',1)
plot(alpha_TAT,Cm_ac_TAT_Fg(4,:),'r--','LineWidth',1)
plot(alpha_TAT,Cm_ac_TAT_F(4,:),'b--','LineWidth',1)
plot(alpha_TAT,Cm_ac_TAT_Fg(7,:),'r-.','LineWidth',1)
plot(alpha_TAT,Cm_ac_TAT_F(7,:),'b-.','LineWidth',1)
hold off
grid on
legend('Gurney \alpha_F=-5\circ','Trad \alpha_F=-5\circ','Gurney \alpha_F=2.5\circ','Trad \alpha_F=2.5\circ','Gurney \alpha_F=7.5\circ','Trad \alpha_F=7.5\circ')
xlabel('\alpha [deg]')
ylabel('c_{m,ac}')
set(gca,'TickLabelInterpreter','latex')
set(gcf,'color','w');
saveas(fig_cmac_flap_gurney, 'fig_cmac_flap_gurney', 'png')

%%
%%%%%%%%%%%%%%%% ERRORS %%%%%%%%%%%%%%%%


%% Error on Cl

% Calculating

% Kutta-Joukowski lift coeff

index_0 = find(rad2deg(alpha_jouk)==0); % index on alpha = 0°
index_8 = find(rad2deg(alpha_jouk)==8); % index on alpha = 8° (region of stall)


cl_exp_interp_JK = interp1(alpha_exp,cl_exp,rad2deg(alpha_jouk)); % interpolation of the experimental data on 

err_KJ = mean(abs(cl_exp_interp_JK(index_0:index_8)-cl_kj(index_0:index_8))./abs(cl_exp_interp_JK(index_0:index_8)))*100 %[%]

% Basic Aerodynamics lift coeff

err_BA = mean(abs(cl_exp_interp_JK(index_0:index_8)-cl_ba(index_0:index_8))./abs(cl_exp_interp_JK(index_0:index_8)))*100 %[%]

% Perfect Mapping lift coeff

index_0_100 = 24; % index on alpha = 0°
index_8_100 = 54; % index on alpha = 8° (region of stall)

err_PM = mean(abs(cl_exp_PFT(index_0_100:index_8_100)-Cl_map(index_0_100:index_8_100))./abs(cl_exp_PFT(index_0_100:index_8_100)))*100 %[%]

% Thin Airfoil Theory

index_0_TAT = find(alpha_TAT == 0); % index on alpha = 0°
index_8_TAT = find(alpha_TAT == 8);

err_TAT = mean(abs(Cl_exp_TAT(index_0_TAT:index_8_TAT)-Cl_TAT(index_0_TAT:index_8_TAT))./abs(Cl_exp_TAT(index_0_TAT:index_8_TAT)))*100 %[%]


%% Error on dCl/d(alpha)

% Experimental slope

p = polyfit(rad2deg(alpha_jouk(index_0:index_8)),cl_exp_interp_JK(index_0:index_8),1);
dcl_exp = p(1);

% Kutta-Joukowski slope

p = polyfit(rad2deg(alpha_jouk(index_0:index_8)),cl_kj(index_0:index_8),1);
dcl_KJ = p(1);

% Basic Aerodynamics slope

p = polyfit(rad2deg(alpha_jouk(index_0:index_8)),cl_ba(index_0:index_8),1);
dcl_BA = p(1);

% Perfect Mapping slope

p = polyfit(rad2deg(alpha3(index_0_100:index_8_100)),Cl_map(index_0_100:index_8_100),1);
dcl_map = p(1);

% Thin Airfoil Theory

p = polyfit(alpha_TAT(index_0_TAT:index_8_TAT),Cl_TAT(index_0_TAT:index_8_TAT),1);
dcl_TAT = p(1);

% Errors

err_dcl_KJ = abs(dcl_exp-dcl_KJ)/dcl_exp*100

err_dcl_BA = abs(dcl_exp-dcl_BA)/dcl_exp*100

err_dcl_map = abs(dcl_exp-dcl_map)/dcl_exp*100

err_dcl_TAT = abs(dcl_exp-dcl_TAT)/dcl_exp*100