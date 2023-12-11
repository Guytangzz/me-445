% Parameters
a = 1; % Parameter for Joukowski transform
r = 0.5; % Radius of the cylinder
theta = linspace(0, 2*pi, 100); % Angles for points around the cylinder

% Complex variable for points on the cylinder
z_cylinder = r * exp(1i * theta);

% Joukowski transform
w = z_cylinder + (a^2) ./ z_cylinder;

% Plotting
figure;

% Plot the cylinder
subplot(1, 2, 1);
plot(real(z_cylinder), imag(z_cylinder), 'b');
title('Cylinder');
xlabel('Real axis');
ylabel('Imaginary axis');
axis equal;

% Plot the airfoil using Joukowski transform
subplot(1, 2, 2);
plot(real(w), imag(w), 'r');
title('Airfoil (Joukowski Transform)');
xlabel('Real axis');
ylabel('Imaginary axis');
axis equal;

% Add a circle to represent the transformed cylinder
hold on;
circle_points = exp(1i * linspace(0, 2*pi, 100));
plot(a * real(circle_points), a * imag(circle_points), '--k');
hold off;
legend('Transformed Airfoil', 'Original Cylinder');

sgtitle('Joukowski Transform: Cylinder to Airfoil');
