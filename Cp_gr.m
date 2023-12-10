function [Cp] = Cp_gr(Uinf,k,alpha,U_gr,alpha_gr)
% this function evaluate the Cp at each point location by taking the ground
% effect and the infinit speed
% The total infinit speed or angle represent the virtual speed that the
% airfoil see when it is close to the ground

    U_inf_gr_h = Uinf.*cos(alpha) + U_gr.*cos(alpha+alpha_gr) ; % horizontal fluid speed seen by the airfoil
    U_inf_gr_v = Uinf.*sin(alpha) + U_gr.*sin(alpha+alpha_gr) ; % vertical fluid speed seen by the airfoil
    delta_u = (1/2).*k ;                                       % local jump in tangiential velocity
    U_inf_gr = sqrt(U_inf_gr_h.^2 + U_inf_gr_v.^2) ;            % total infinit speed with ground effect
    alpha_inf_gr = acos(U_inf_gr_h/U_inf_gr) ;                  % total infinit angle with ground effect
    U_up = sqrt((U_inf_gr.*cos(alpha_inf_gr)+delta_u).^2 + (U_inf_gr.*sin(alpha_inf_gr)).^2) ; % upper speed
    U_low = sqrt((U_inf_gr.*cos(alpha_inf_gr)-delta_u).^2 + (U_inf_gr.*sin(alpha_inf_gr)).^2) ; % lower speed
    Cpu = 1 - (U_up.^2 ./ U_inf_gr.^2) ;                        % upper pression coef with ground effect
    Cpl = 1 - (U_low.^2 ./ U_inf_gr.^2) ;                       % lower pression coef with ground effect
    Cp = Cpl - Cpu ;                                            % pression coef with ground effect

end
