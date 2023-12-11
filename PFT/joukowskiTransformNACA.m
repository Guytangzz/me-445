function [z_cylinder, w_cylinder, z_airfoil, w_airfoil] = joukowskiTransformNACA(naca)

    % Parse NACA 4-digit parameters
    a = str2double(naca(1))/100;    % Maximum camber percentage
    b = str2double(naca(2))/10;     % Location of maximum camber
    c = str2double(naca(3:4))/100;  % Thickness percentage

    % Define cylinder in z-plane
    theta_cylinder = linspace(0, 2*pi, 1000);
    z_cylinder = exp(1i * theta_cylinder);

    % Joukowski transform for cylinder
    w_cylinder = z_cylinder + (a^2) ./ z_cylinder;

    % Calculate NACA 4-digit camber line
    theta_airfoil = linspace(0, 2*pi, 1000);
    x_camber = 0.5 * (1 - cos(theta_airfoil));
    y_camber = a * x_camber .* (1 - x_camber/b) ./ (b^2);

    % Thickness distribution
    yt = 5 * c * (0.2969 * sqrt(x_camber) - 0.1260 * x_camber - 0.3516 * x_camber.^2 + 0.2843 * x_camber.^3 - 0.1015 * x_camber.^4);

    % Apply Joukowski transform to NACA airfoil
    z_airfoil = x_camber + 1i * (y_camber + yt);
    w_airfoil = z_airfoil + (a^2) ./ z_airfoil;

    % Plotting
    figure;

    % Plot the original cylinder in z-plane
    subplot(1, 2, 1);
    plot(real(z_cylinder), imag(z_cylinder), 'b');
    title('Original Cylinder (z-plane)');
    xlabel('Real axis');
    ylabel('Imaginary axis');
    axis equal;

    % Plot the transformed NACA airfoil in the real plane
    subplot(1, 2, 2);
    plot(real(w_airfoil), imag(w_airfoil), 'r');
    title('Transformed NACA Airfoil (Real plane)');
    xlabel('Real axis');
    ylabel('Imaginary axis');
    axis equal;

    sgtitle(['Joukowski Transform: NACA ' naca ' Airfoil']);
end
