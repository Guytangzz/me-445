% Sinclair Augereau --/11/2023
% ground effect approch with thin airfoil theory
clear all
close all

%% definition of the variables
nb_alpha = 100 ;                      % [-], number of angle alpha evaluated
alpha_deg = linspace(-6,22,nb_alpha); % angle of attack in degree
alpha_rad = alpha_deg *pi/180 ;

% NACA4412 profile : but you can choose any other 4-digit of your choice
m = 0.04;                             % [-], relative maximum camber
p = 0.4;                              % [-], relative location of maximum camber
t = 0.12;                             % [-], relative thikness
c = 1;                                % [m], cord length
Uinf = 10 ;                           % [m/s], speed of the airflow in amout

CL = lift_coef(m, p, c, alpha_rad); % hand computation
% [A0, An] = A_n_computation(m, p, 1) ; % code computation
% CL_code = 2.*pi.*alpha_rad + pi.*(An(1)-2*A0); 
% % to check if our computation of A0 and A1 was good 

%% get the point of the naca profile selcted by using this function
nb = 100 ;          % [-], number of point of the discrétisation of the corde
alpha = 10*pi/180;  % [rad], Angle of attack (AOA)
h = 0.5;            % [m], height with respect to ground (mesured at the trailing edge)

% Naca discretised points with : 
[Z , Zl , Zu, k] = naca_point(m, p, t, c, alpha, h, nb, Uinf) ;
                   % Z the cord points
                   % Zl the lower profile points
                   % Zu the upper profile points
                   % k dicretised vortex strength

% definition of the symetrical point of the cord
Zsy = [Z(1,:) ; - Z(2,:)] ;


%% computation of the Cl without the ground but with the discretized methode
Cl_inf = zeros(nb_alpha,1) ;                                                                % initalisation
for i =1:nb_alpha                                                                                                                                       
    [aa , bb , cc, k_inf] = naca_point(m, p, t, c, alpha_rad(i), h, nb, Uinf) ;                % NACA discretised
    Cpu_inf = - abs(k_inf)/Uinf ;                                                           % [-], upper pression coef
    Cpl_inf = abs(k_inf)/Uinf ;                                                             % [-], lower pression coef
    Cl_inf(i) = sum(Cpl_inf(1,2:length(Cpl_inf)) - Cpu_inf(1,2:length(Cpu_inf)))*(1/nb) ;   % [-], Lift coef
    % Cl_test(i) = 1/Uinf * integral ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on
plot(alpha_deg, CL);
plot(alpha_deg,Cl_inf,"-",'Color','red'); % plot the Cl values found with the discretized methode
% plot(alpha_deg, CL_code,"-","Color","green");
title(sprintf('NACA %d%d%d \n Lift coeficient without ground effect', m*100, p*10, t*100));
legend('Thin airfoil Theory',sprintf('Thin airfoil Theory discretised in %d points',nb),'FontSize', 14);
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
title(sprintf('NACA %d%d%d \n AOA = %d°  and  h/c = %.2f', m*100, p*10, t*100, alpha*180/pi, h/c), 'FontSize', 14);
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

