% Sinclair Augereau --/11/2023
% ground effect approch with thin airfoil theory
clear all
close all

% definition of the variables
alpha_deg = linspace(-6,22,100); % angle of attack in degree
alpha_rad = alpha_deg *pi/180 ;
% NACA4412 profile : but you can choose any other 4-digit of your choice
m = 0.04; % relative maximum camber
p = 0.4; % relative location of maximum camber
t = 0.12; % relative thikness
c = 1;
Uinf = 10 ; %[m/s], speed of the airflow in amout

CL = lift_coef(m, p, c, alpha_rad);

%% get the point of the naca profile selcted by using this function
nb = 100 ; % number of point of the discrétisation of the corde
alpha = 10*pi/180; % Angle of attack in rad (AOA)
h = 0.5; % high with respect to ground (mesure at the trelling edge)
[Z , Zl , Zu, k] = naca_point(m, p, t, c, alpha, h, nb, Uinf) ; % points with Z the corde point; Zl the lower profile point and Zu the upper profile point

% definition of the symetrical point of the cord
Zsy = [Z(1,:) ; - Z(2,:)] ;


%% computation of the Cl without the ground but with the discretized methode
for i =1:20 % variation of the AOA from 1 to 20°
    alpha_t(i) = i*pi/180; % Angle of attack in rad (AOA)
    [aa , bb , cc, k_inf] = naca_point(m, p, t, c, alpha_t(i), h, nb, Uinf) ;
    Cpu_inf = - abs(k_inf)/Uinf ;
    Cpl_inf = abs(k_inf)/Uinf ;
    Cl_inf(i) = sum(Cpl_inf(1,2:length(Cpl_inf)) - Cpu_inf(1,2:length(Cpu_inf)))*(1/nb) ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on
plot(alpha_deg, CL);
plot(alpha_t*180/pi,Cl_inf,"d",'Color','red'); % plot the Cl values found with the discretized methode
title('Lift coeficient of a airfoil without ground effect');
xlabel('alpha [°]');
ylabel('CL [-]');
hold off

%% graphic of the NACA 4412 at a high of h/c
figure;
hold on
plot(Z(1,:),Z(2,:),'.') % cord
plot(Zl(1,:),Zl(2,:)) % lower
plot(Zu(1,:),Zu(2,:)) % upper
for i=1:length(Zl(1,:))
    plot([Zu(1,i), Zl(1,i)] , [Zu(2,i), Zl(2,i)])
end
% plot the symetrical corde line :
plot(Zsy(1,:),Zsy(2,:),'.') % symetrical cord
% plot a line to schematized the ground :
plot([0, Z(1,end)] , [0, 0])
%title('NACA 4412 profile discretized');
title(sprintf('NACA %d%d%d \n AOA = %d°  and  h/c = %d', m*100, p*10, t*100, alpha*180/pi, h/c), 'FontSize', 14);
xlabel('x/c [-]','FontSize', 14);
ylabel('y/c [-]','FontSize', 14);
axis equal;
hold off

%% graphic of the vortex strength dicretization on the corde with respect to the naca profile
figure;
%plot(Z(1,(1:length(Z)-1)), k((1:length(Z)-1))) % we don't plot the last point as it is a singularity
plot(Z(1,:), k)
title(sprintf('NACA %d%d%d', m*100, p*10, t*100), 'FontSize', 14);
xlabel('x/c [-]','FontSize', 14);
ylabel('vortex strength [circulation/length]','FontSize', 14);

