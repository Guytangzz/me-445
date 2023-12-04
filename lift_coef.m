function CL = lift_coef(m, p, c, alpha)
% computation of the lift coefition of a naca profile
% by using the thine arifoil theory
% NACAmpt 
% alpha is supposed to be a vector which represent the angle of atack


   %  % range of validity
   %  % p/c = (1-cos(theta))/2  <=> theta = acos(1-2*p/c)
   %  theta_mid = acos(1-2*p/c);
     
   % % front = (theta >= 0) & (x < theta_mid);
   % % back = (theta >= theta_mid) & (theta <= pi);
   % 
   %  % for front case
   %  yc_dot_front1 = @(theta)  (2*m/(p^2)).*(p-(c.*(1-cos(theta))/2));
   %  yc_dot_front0 = @(theta)  (2*m/(p^2)).*(p-(c.*(1-cos(theta))/2)).*cos(theta);
   % 
   %  % for back case
   %  yc_dot_back1 = @(theta)  (2*m/((1-p)^2)).*(p-(c*(1-cos(theta))/2));
   %  yc_dot_back0 = @(theta)  (2*m/((1-p)^2)).*(p-(c*(1-cos(theta))/2)).*cos(theta);
   % 
   %  % computation of An
   %  A1_front = (2/pi) .* integral(yc_dot_front1, 0,theta_mid);
   %  A0_front = (2/pi) .* integral(yc_dot_front0, 0, theta_mid);
   %  A1_back = (2/pi) .* integral(yc_dot_back1,theta_mid, pi);
   %  A0_back = (2/pi) .* integral(yc_dot_back0, theta_mid, pi);
   % 
   %  A0 = A0_front + A0_back ;
   %  A1 = A1_front + A1_back ;

    % MANUAL COMPUTATION
    theta_p = acos(1-2*p/c);
    A0 = (m/(pi*p^2)) * ((2*p-1)*theta_p + sin(theta_p)) + (m/(pi*(1-p)^2))*((2*p-1)*(pi-theta_p)-sin(theta_p)) ;
    A1 = (2*m/(pi*p^2)) * ((2*p-1)*sin(theta_p) + 1/4 *sin(2*theta_p) + theta_p/2) - (2*m/(pi*(1-p)^2)) * ((2*p-1)*sin(theta_p) + 1/4 *sin(2*theta_p) - 1/2 *(pi-theta_p)) ;
   
    % compitation of CL
    CL = 2.*pi.*alpha + pi.*(A1-2*A0);
end