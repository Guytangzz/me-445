function k = vortex(xc_front, xc_back, m, p, alpha, Uinf, theta_front, theta_back)
% this programe aims to evalutae the vorteces strenght at each location
% specified of a discretized naca profile

% Transforme the xc coordinate into angle of a cercle (with 0 at the
% ledding edge and pi at the trelling edge
    theta_xc_front = acos(1-2.*xc_front);
    theta_xc_back = acos(1-2.*xc_back);

% computation of the effective angle pf attack (angle between the cord and
% the incident horizontal flow
    alpha_eff_front = -(theta_front-alpha);
    alpha_eff_back = -(theta_back-alpha);

% Compute A0 and An
    n = 25 ; % number of term in the Fourier transform
    [A0, An] = A_n_computation(m, p, n) ;

% computation of k1 :
    k1_front = 2.*Uinf.*(alpha_eff_front-A0).*(1./tan(theta_xc_front)) ;
    k1_back = 2.*Uinf.*(alpha_eff_back-A0).*(1./tan(theta_xc_back)) ;

% computation of k2 :
    k2_front = 2 * Uinf .* (alpha_eff_front - A0)./sin(theta_xc_front) ;
    k2_back = 2 * Uinf .* (alpha_eff_back - A0)./sin(theta_xc_back) ;

% computation of k3 : 
    k3_front = 2.*Uinf.* sin(theta_xc_front)*sum(An); 
    k3_back = 2.*Uinf.* sin(theta_xc_back)*sum(An);
    %k3_front = 2.*Uinf.* (m*sin(theta_xc_front))./(p^2) ; 
    %k3_back = 2.*Uinf.* (m*sin(theta_xc_back))./((1-p)^2) ;

% get back the total k :
    k_front = k1_front + k2_front + k3_front ;
    k_back = k1_back + k2_back + k3_back ;

    k =  [k_front, k_back] ;
end