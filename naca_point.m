function [Z , Zl , Zu, k] = naca_point(m, p, t, c, alpha, h, nb, Uinf)
% Sinclair Augereau --/11/2023
% function that define all point of a naca profile
% it also uses the vortex function to compute the vorteces strenght at each location
% specified of the discretized naca profile
% alpha is a unique value here in radian
% NACAmpT ex : NACA4412 m = 0.04 ; p = 0.4 ; t = 0.12
% h is the high of the airfoil (distance from the treling edge to the ground) 
% nb is the number of point discretisation

%% all is adimentionalized with respect to the corde c

% repartition of the point evenly between front and back
    nb_front = round(nb*p) ; % we take the round to get a integer
    nb_back = nb - nb_front ;

% discretisation of the horizontal position of the point with a alpha=0
% separated with front and back point
    xc_front = linspace(0,p,nb_front); 
    xc_back = linspace(p+p/nb_front,c,nb_back); 

% computation of the high of the cordline points
% this equation is a naca standardization
    yc_front = (m/(p^2)) .* (2.*p.*xc_front - xc_front.^2);
    yc_back = (m/((1-p)^2)) .* ((1-2*p) + 2.*p.*xc_back - xc_back.^2);

%     % debug
%     figure;
%     plot(xc_front,yc_front),hold on
%     plot(xc_back,yc_back),hold off
% 
%     for i = 2:nb_front
%         dyc_front_test(i) = (yc_front(i) - yc_front(i-1)) ./ (xc_front(i) - xc_front(i-1));
%     end
% 
%     for i = 2:nb_back
%         dyc_back_test(i) = (yc_back(i) - yc_back(i-1)) ./ (xc_back(i) - xc_back(i-1));
%     end
% 
%     figure;
%     plot(xc_front,dyc_front_test), hold on 
%     plot(xc_back,dyc_back_test), hold off



% know we take the angle of attack into account
    x0_front = 1 - (1-xc_front)*cos(alpha) ;
    x0_back = 1 - (1-xc_back)*cos(alpha) ;
    y0_front = yc_front + (1-xc_front)*sin(alpha) + h/c ;
    y0_back = yc_back + (1-xc_back)*sin(alpha) + h/c ;

    Z = [ [x0_front, x0_back] ; [y0_front, y0_back]] ;

%% know here is computed the lower and upper point of the profile

% first the thinkness is computed :
% constant of naca standardization
    a0 = 0.2969 ;
    a1 = -0.126 ;
    a2 = -0.3516 ;
    a3 = 0.2843 ;
    a4 = -0.1015 ; % or -0.1036 for  a closed trailing edge

    yt_front = t/0.2 * (a0.*xc_front.^0.5 + a1.*xc_front.^1 + a2.*xc_front.^2 + a3.*xc_front.^3 + a4.*xc_front.^4);
    yt_back = t/0.2 * (a0.*xc_back.^0.5 + a1.*xc_back.^1 + a2.*xc_back.^2 + a3.*xc_back.^3 + a4.*xc_back.^4);

% Then the derivative of the profile :
    dyc_front = (2*m/(p^2)) .* (p-xc_front);
    dyc_back = (2*m/((1-p)^2)) .* (p-xc_back);

% definition of a useful angles
    theta_front = atan(dyc_front) ;
    theta_back = atan(dyc_back) ;
    xi_front = theta_front - alpha;
    xi_back = theta_back - alpha;

%% know here is computed the lower point of the profile   
% definition of the x lower point of the profile
    xl_front = x0_front + yt_front.*sin(xi_front);
    xl_back = x0_back + yt_back.*sin(xi_back);

% definition of the y lower point of the profile
    yl_front = y0_front - yt_front.*cos(xi_front);
    yl_back = y0_back - yt_back.*cos(xi_back);

% get the values back
    Zl = [ [xl_front, xl_back] ; [yl_front, yl_back]] ;


%% know here is computed the upper point of the profile
% definition of the x lower point of the profile
    xu_front = x0_front - yt_front.*sin(xi_front);
    xu_back = x0_back - yt_back.*sin(xi_back);

% definition of the y lower point of the profile
    yu_front = y0_front + yt_front.*cos(xi_front);
    yu_back = y0_back + yt_back.*cos(xi_back);

% get the values back
    Zu = [ [xu_front, xu_back] ; [yu_front, yu_back]] ;

% compute the vorteces strenght at each location
% specified of a discretized naca profile
    k = vortex(xc_front, xc_back, m, p, alpha, Uinf, theta_front, theta_back);

    
    % for i=1:length(theta_front)
    %     Cl_front = Cl_front + (1/Uinf * k*sin(theta_front)*(theta_front(i+1)-theta_front(i)) );
    % end
    %     Cl_back
    % Cl_test = 1/Uinf * k*sin(,0,pi) ;
end