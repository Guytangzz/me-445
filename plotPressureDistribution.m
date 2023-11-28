function plotPressureDistribution(airfoilFile, chordLength)
    % Step 1: Read the airfoil geometry
    airfoilData = dlmread(airfoilFile, ' ', 1, 0);

    % Extract x and y coordinates
    x = airfoilData(:, 1);
    y = airfoilData(:, 2);

    % Step 2: Apply the Joukowsky Transformation
    z = x + 1i * y;
    w = z + (chordLength^2) ./ z;

    % Step 3: Calculate the Pressure Distribution
    dw_dz = 1 - (chordLength^2) ./ (z.^2);
    Cp = 1 - abs(dw_dz).^2;

    % Step 4: Plot the Pressure Distribution for both sides
    figure;
    
    % Lower surface (negative y-values)
    plot(x, Cp, 'b-', 'LineWidth', 2);
    hold on;

    % Upper surface (positive y-values)
    plot(x, Cp, 'r-', 'LineWidth', 2);

    title('Pressure Distribution over Airfoil');
    xlabel('x');
    ylabel('Cp');
    legend('Lower Surface', 'Upper Surface');
    grid on;
    hold off;
end
