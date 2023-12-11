% Ground Effect Modeling using Potential Flow Theory with Panel Method

% Parameters
U_inf = 1.0;  % Freestream velocity
alpha = 5;    % Angle of attack in degrees
h = 0.1;      % Ground clearance

% Airfoil coordinates (example: NACA 0012)
[x_airfoil, y_airfoil] = naca_airfoil('0012', 100);

% Extend airfoil to create a closed body
x_body = [x_airfoil; flipud(x_airfoil)];
y_body = [y_airfoil; -flipud(y_airfoil)];

% Plot the body
figure;
plot(x_body, y_body, 'k-', 'LineWidth', 2);
hold on;
plot(x_airfoil, y_airfoil, 'r-', 'LineWidth', 1);
xlabel('X');
ylabel('Y');
title('Airfoil with Ground Effect');

% Calculate flow parameters
alpha_rad = deg2rad(alpha);
V_inf_x = U_inf * cos(alpha_rad);
V_inf_y = U_inf * sin(alpha_rad);

% Calculate image source strength for ground effect
gamma_ground = -2 * pi * h * U_inf;

% Add image source below the airfoil for ground effect
x_image = x_airfoil;
y_image = -y_airfoil - h;
gamma_image = gamma_ground * ones(size(x_image));

% Combine actual and image bodies
x_combined = [x_airfoil; x_image];
y_combined = [y_airfoil; y_image];
gamma_combined = [zeros(size(x_airfoil)); gamma_image];

% Solve for velocity potential using panel method
[phi, dphi_dx, dphi_dy] = panel_method(x_combined, y_combined, gamma_combined, V_inf_x, V_inf_y);

% Calculate pressure coefficient
cp = 1 - (sqrt(dphi_dx.^2 + dphi_dy.^2)) / U_inf;

% Calculate lift and drag
length_panels = length(x_airfoil);
lift = sum(cp(1:length_panels-1) .* diff(x_airfoil));
drag = sum(cp(1:length_panels-1) .* diff(y_airfoil));

fprintf('Lift: %.4f, Drag: %.4f\n', lift, drag);

fprintf('Lift: %.4f, Drag: %.4f\n', lift, drag);

% Plot pressure coefficient distribution
figure;
plot(x_airfoil, cp(1:length_panels), 'b-', 'LineWidth', 2);
xlabel('X');
ylabel('C_p');
title('Pressure Coefficient Distribution');

% Additional functions

function [x, y] = naca_airfoil(naca, num_points)
    % Generate NACA airfoil coordinates
    % Example: naca_airfoil('0012', 100)
    % Returns x and y coordinates for a NACA 0012 airfoil with 100 points.

    t = str2double(naca(3)) / 100;  % Maximum thickness
    m = str2double(naca(1)) / 100;  % Maximum camber
    p = str2double(naca(2)) / 10;   % Location of maximum camber

    x = linspace(0, 1, num_points);
    y_camber = m * (x - p).^2 / (p * (1 - p)^2);
    y_thickness = 5 * t * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x.^2 + 0.2843 * x.^3 - 0.1015 * x.^4);

    y_upper = y_camber + 0.5 * y_thickness;
    y_lower = y_camber - 0.5 * y_thickness;

    x = [x, fliplr(x)];
    y = [y_upper, fliplr(y_lower)];
end

function [phi, dphi_dx, dphi_dy] = panel_method(x, y, gamma, V_inf_x, V_inf_y)
    % Panel method for potential flow analysis
    % Returns the velocity potential (phi) and its derivatives (dphi_dx, dphi_dy).

    num_panels = length(x) - 1;

    % Initialize matrices for linear system
    A = zeros(num_panels, num_panels);
    B = zeros(num_panels, 1);

    % Loop through each panel
    for i = 1:num_panels
        xi = x(i);
        yi = y(i);
        xi1 = x(i + 1);
        yi1 = y(i + 1);

        % Panel length and orientation
        ds = sqrt((xi1 - xi)^2 + (yi1 - yi)^2);
        nx = -(yi1 - yi) / ds;  % Outward normal x-component
        ny = (xi1 - xi) / ds;   % Outward normal y-component

        % Influence coefficients from panel j on panel i
        for j = 1:num_panels
            xj = x(j);
            yj = y(j);
            xj1 = x(j + 1);
            yj1 = y(j + 1);

            % Panel center
            xc = 0.5 * (xi + xi1);
            yc = 0.5 * (yi + yi1);

            % Panel coordinates
            xij = xc - xj;
            yij = yc - yj;
            xij1 = xc - xj1;
            yij1 = yc - yj1;

            % Panel normal
            njx = -(yj1 - yj) / ds;
            njy = (xj1 - xj) / ds;

            % Compute influence coefficients
            A(i, j) = -2 * pi * (ny * (atan2(yij1, xij1) - atan2(yij1, xij)) - ny * (atan2(yij, xij1) - atan2(yij, xij)));
        end

        % Right-hand side of the linear system
        B(i) = -V_inf_x * nx - V_inf_y * ny;

        % Contribution from the vortex sheet
        for j = 1:num_panels
            B(i) = B(i) - gamma(j) / (2 * pi) * (ny * (atan2(yij1, xij1) - atan2(yij, xij1)) - ny * (atan2(yij1, xij) - atan2(yij, xij)));
        end
    end

    % Solve the linear system to find panel strengths (gamma)
    gamma = linsolve(A, B);

    % Compute velocity potential and its derivatives
    phi = zeros(size(x));
    dphi_dx = zeros(size(x));
    dphi_dy = zeros(size(x));

    for i = 1:num_panels
        for j = 1:num_panels
            xj = x(j);
            yj = y(j);
            xj1 = x(j + 1);
            yj1 = y(j + 1);

            % Panel center
            xc = 0.5 * (xi + xi1);
            yc = 0.5 * (yi + yi1);

            % Panel coordinates
            xij = xc - xj;
            yij = yc - yj;
            xij1 = xc - xj1;
            yij1 = yc - yj1;

            % Panel normal
            njx = -(yj1 - yj) / ds;
            njy = (xj1 - xj) / ds;

            % Velocity potential and derivatives due to panel j
            phi(i) = phi(i) + gamma(j) / (2 * pi) * (ny * (atan2(yij1, xij1) - atan2(yij, xij)) - ny * (atan2(yij1, xij) - atan2(yij, xij)));
            dphi_dx(i) = dphi_dx(i) + gamma(j) / (2 * pi) * (njx * (log((xij1)^2 + (yij1)^2) - log((xij)^2 + (yij)^2)));
            dphi_dy(i) = dphi_dy(i) + gamma(j) / (2 * pi) * (njy * (log((xij1)^2 + (yij1)^2) - log((xij)^2 + (yij)^2)));
        end
    end

    % Superimpose freestream contribution
    phi = phi + V_inf_x * x + V_inf_y * y;
    dphi_dx = dphi_dx + V_inf_x;
    dphi_dy = dphi_dy + V_inf_y;
end
