function [CL] = Cp_gr_curve(m, p, t, c, alpha, h, nb, Uinf)
% here the Cl is compute for varius AOA and it also take the ground effect

dxc = 1/nb ;        % [-], adimentionalized infinetesimal discretized length
% initialisation
    CL = zeros(1,length(alpha));
    for i=1:length(alpha)
        % Naca discretised points with :
        [Z , Zl , Zu, k] = naca_point(m, p, t, c, alpha(i), h, nb, Uinf) ;
        % compute the strength and direction of the induce speed of a vortex over a point 
        [U_gr , alpha_gr] = vortex_induce_speed (Z,k,dxc) ; % call of the compute function
        % compute the Cp with ground effect
        [Cp] = Cp_gr(Uinf,k,alpha(i),U_gr,alpha_gr) ; % call of the compute function
        % compute the Cl induce by a Cp
        CL(i) = sum(Cp.*dxc) ; % discretized integral
    end

end