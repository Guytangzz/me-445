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

% computation of A0 :
% A0 = (1/pi) * integral from 0 to pi of (dyc/dx) dtheta witch is
% equivalent to (for a naca case) :
    A0 = (2*m/(pi)) * ( (1/(2*p^2)) * ( 2*sqrt(-(p-1)) + (2*p-1)*acos(1-2*p) ) + ( (1/(2*(1-p)^2)) * (pi*(2*p-1) - 2*sqrt(-(p-1)) + (1-2*p)*acos(1-2*p)) ) ) ;

% computation of k1 :
    k1_front = 2.*Uinf.*(alpha_eff_front-A0).*(1./tan(theta_xc_front)) ;
    k1_back = 2.*Uinf.*(alpha_eff_back-A0).*(1./tan(theta_xc_back)) ;

% computation of k2 :
    k2_front = 2 * Uinf .* (alpha_eff_front - A0) ;
    k2_back = 2 * Uinf .* (alpha_eff_back - A0) ;

% computation of k3 :  
    k3_front = 2.*Uinf.* (m*sin(theta_xc_front))./(p^2) ; 
    k3_back = 2.*Uinf.* (m*sin(theta_xc_back))./((1-p)^2) ;

% get back the total k :
    k_front = k1_front + k2_front ; %+ k3_front ;
    k_back = k1_back + k2_back ; %+ k3_back ;

    k =  [k_front, k_back] ;
end