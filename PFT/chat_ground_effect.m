%% Ground Effect Simulation
Gamma_ground_effect = -2 * pi * U_inf * c;  % Circulation to simulate ground effect

W_ground_effect = (U_inf * exp(-1i * alpha_jouk) + 1i * Gamma_joukowski / (2 * pi * (zeta_circ - zeta_0)) - U_inf * R_joukowski^2 * exp(1i * alpha_jouk) / ((zeta_circ - zeta_0).^2) + 1i * Gamma_ground_effect / (2 * pi * (zeta_circ - zeta_0)))./(1 - a^2 ./ (zeta_circ.^2));

u_ground_effect = real(W_ground_effect);
v_ground_effect = -imag(W_ground_effect);
V_ground_effect = sqrt(u_ground_effect.^2 + v_ground_effect.^2);
Cp_ground_effect = 1 - V_ground_effect.^2 / U_inf^2;

% Adjust dimensions to match length(theta)
Cp_ground_effect_upper = flip(Cp_ground_effect(1:length(theta)/2));
Cp_ground_effect_lower = Cp_ground_effect(length(theta)/2+1:end);

c_n_ba(i) = trapz(x_down_approx, Cp_ground_effect_lower) - trapz(x_up_approx, Cp_ground_effect_upper);
cl_ba(i) = cos(alpha_jouk(i)) * c_n_ba(i);

L_kj = rho * U_inf * (Gamma_joukowski(i) + Gamma_ground_effect);
cl_kj(i) = 2 * L_kj / (rho * U_inf^2 * c);

cm_le_ba(i) = trapz(x_up_approx, Cp_ground_effect_upper .* x_up_approx') - trapz(x_down_approx, Cp_ground_effect_lower .* x_down_approx');
cm_14_ground_effect_i = cm_le_ba(i) + 0.25 * c_n_ba(i);



% Plot Cp distribution with ground effect
figure;
plot(real(z/c), Cp_ground_effect, '-b', 'LineWidth', 2);
xlabel('Chord-wise Coordinate (x/c)');
ylabel('Pressure Coefficient (Cp)');
title('Joukowski Airfoil Cp Distribution with Ground Effect');
grid on;

%% Coefficients computation with ground effect

% Coefficients computation with ground effect
for i = 1:length(alpha_jouk)
    W_ground_effect_i = (U_inf * exp(-1i * alpha_jouk(i)) + 1i * Gamma_joukowski(i) / (2 * pi * (zeta_circ - zeta_0)) - U_inf * R_joukowski^2 * exp(1i * alpha_jouk(i)) / ((zeta_circ - zeta_0).^2) + 1i * Gamma_ground_effect / (2 * pi * (zeta_circ - zeta_0)))./(1 - a^2 ./ (zeta_circ.^2));
    u_ground_effect_i = real(W_ground_effect_i);
    v_ground_effect_i = -imag(W_ground_effect_i);
    V_ground_effect_i = sqrt(u_ground_effect_i.^2 + v_ground_effect_i.^2);
    Cp_ground_effect_i = 1 - V_ground_effect_i.^2 / U_inf^2;

    Cp_up_ground_effect_i = flip(Cp_ground_effect_i(1:length(theta)/2));
    Cp_down_ground_effect_i = Cp_ground_effect_i(length(theta)/2+1:end);

    c_n_ba(i) = trapz(x_down_approx, Cp_down_ground_effect_i) - trapz(x_up_approx, Cp_up_ground_effect_i);
    cl_ba(i) = cos(alpha_jouk(i)) * c_n_ba(i);

    L_kj = rho * U_inf * (Gamma_joukowski(i) + Gamma_ground_effect);
    cl_kj(i) = 2 * L_kj / (rho * U_inf^2 * c);

    cm_le_ba(i) = trapz(x_up_approx, Cp_up_ground_effect_i .* x_up_approx') - trapz(x_down_approx, Cp_down_ground_effect_i .* x_down_approx');
    cm_14_ground_effect_i = cm_le_ba(i) + 0.25 * c_n_ba(i);
end

%% RESULTS WITH GROUND EFFECT
