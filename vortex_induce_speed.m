function [U,alpha] = vortex_induce_speed (Z,k,dxc)
% This function compute for each cord point the strength and direction of 
% the induce speed due to all the symetrical vortex of the cord line 

% Uh is a vecteur that give the horizontal induce speed for each point
% Uv is a vecteur that give the vertical induce speed for each point
    
% definition of the symetrical point of the cord
    Zsy = [Z(1,:) ; - Z(2,:)] ;

% initialisation
Uh_pt = zeros(1,length(Zsy)) ;
Uv_pt = zeros(1,length(Zsy)) ;
Uh = zeros(1,length(Z)) ;
Uv = zeros(1,length(Z)) ;
U = zeros(1,length(Z)) ;
alpha = zeros(1,length(Z)) ;

    for i = 1:length(Z)
        for j = 1:length(Zsy)
            % compute the distance between the point of interest and the center of the vortex of interest
            r = sqrt( (Z(1,i)-Zsy(1,j))^2 + (Z(2,i)-Zsy(2,j))^2)  ; % [m]
            % compute the induce speed
            U_pt = -k(j)*dxc/(2*pi*r) ; % [m/s] we take the negative value of k as vortex of the symetry are of oposite sign of the original airfoil
            alpha_pt = acos((Z(2,i)- Zsy(2,j))/r) ; % [rad]
            Uh_pt(j) = U_pt*cos(alpha_pt) ; % induce horizontal speed
            Uv_pt(j) = U_pt*sin(alpha_pt) ; % induce vertical speed
        end
        Uh(i) = sum(Uh_pt) ; % totale induce horizontal speed on point i due to all symetrical vortex
        Uv(i) = sum(Uv_pt) ; % totale induce vertical speed on point i due to all symetrical vortex
        U(i) = sqrt(Uh(i)^2 + Uv(i)^2) ; % totale induce speed on point i due to all symetrical vortex
        alpha(i) = acos(Uh(i)/U(i)) ; % angle of direction of the total induce speed
    end

end