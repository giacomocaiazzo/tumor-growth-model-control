%% Phase portrait of the reduced error dynamics
clear all
close all
clc

R2 = 0.6;
A12 = 1.0;
A21 = 1.5;

[e1, e2] = meshgrid(-0.4:0.04:0.4, -0.4:0.04:0.4);

de1 = -e1 .* (e1 + e2);
de2 = -(1 + e2) .* (R2 .* e2 + A21 .* e1);

% Normalize the vector field to show directions clearly.
L = sqrt(de1.^2 + de2.^2);
u = de1 ./ (L + eps);
v = de2 ./ (L + eps);

figure('Color', 'w')
hold on
grid on
box on
xlabel('e_1')
ylabel('e_2')
title('Error Dynamics Phase Portrait')
axis([-0.4 0.4 -0.4 0.4])

quiver(e1, e2, u, v, 0.5, 'Color', [0.7 0.7 0.7])

% Sample initial conditions on a circle and an interior grid.
theta = linspace(0, 2*pi, 20);
radius = 0.35;
x_circle = radius * cos(theta);
y_circle = radius * sin(theta);

[xg, yg] = meshgrid(-0.3:0.15:0.3, -0.3:0.15:0.3);
x0_list = [x_circle' y_circle'; xg(:) yg(:)];

for k = 1:size(x0_list, 1)
    sys = @(t,y) [-y(1)*(y(1) + y(2)); 
                  -(1 + y(2))*(R2*y(2) + A21*y(1))];
    
    % Integrate in both time directions to display the local trajectories.
    [~, traj_fwd] = ode45(sys, [0 10], x0_list(k,:));
    [~, traj_bwd] = ode45(sys, [0 -5], x0_list(k,:)); 
    
    plot(traj_fwd(:,1), traj_fwd(:,2), 'b-', 'LineWidth', 1.0)
    plot(traj_bwd(:,1), traj_bwd(:,2), 'b-', 'LineWidth', 1.0)
    
    % Mark the initial point.
    plot(x0_list(k,1), x0_list(k,2), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k')
end

plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'LineWidth', 1.5)
