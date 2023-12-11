function [u, v] = complexPotential(u_inf, alpha, R, Gamma, z, a, zeta_0)
    u = zeros(size(z));
    v = zeros(size(z));
    
    for i = 1:size(z, 1)
        for j = 1:size(z, 2)
            min_ = min(roots([1, -z(i, j), a^2]));
            max_ = max(roots([1, -z(i, j), a^2]));
            if abs(min_ - zeta_0) <= a
                zeta = max_;
            else
                zeta = min_;
            end
            % W_joukowski = (u_inf*exp(-1i*alpha)*...
            %     (zeta - (R.^2*exp(2i*alpha))./(zeta - zeta_0).^2)...
            %     - (1i*Gamma)./(2*pi)./(zeta - zeta_0))...
            %     ./(1-(a.^2)/(zeta.^2));

            W_joukowski = (u_inf.*exp(-1i.*alpha) + 1i.*Gamma./(2*pi.*(zeta-zeta_0)) - u_inf.*R^2.*exp(1i.*alpha)./((zeta-zeta_0).^2))./(1-a^2./(zeta.^2));

            u(i, j) = real(W_joukowski);
            v(i, j) = -imag(W_joukowski);
        end
    end
end
