function [u, v] = complexConjPotential(u_inf, alpha, R, Gamma, z, a, zeta_0, h)
    u = zeros(size(z));
    v = zeros(size(z));
    
    for i = 1:size(z, 1)
        for j = 1:size(z, 2)
            zeta = min(roots([1, -z(i, j), a^2]));
            W_joukowski = conj((u_inf.*exp(-1i.*alpha) - ...
                1i.*Gamma./(2*pi.*(zeta-zeta_0)) - ...
                u_inf.*R^2.*exp(1i.*alpha)./((zeta-zeta_0).^2))...
                ./(1-a^2./((zeta-zeta_0).^2)))...
                -h*i;

            u(i, j) = real(W_joukowski);
            v(i, j) = -imag(W_joukowski);
        end
    end
end