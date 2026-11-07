%% Three-dimensional open-loop phase portrait and equilibria
clear; clc; close all;

R3 = 2.6; R2 = 0.6; D3 = 0.5; k3 = 1.0;
A12 = 1.0; A13 = 2.5; A21 = 1.5; A31 = 0.2;

figure('Color', 'w'); hold on; grid on; box on;
view(135, 25);
xlabel('x_1: Tumor cells'); ylabel('x_2: Hunting cells'); zlabel('x_3: Resting cells');

h1 = plot3(nan, nan, nan, 'Color', [0, 0.25, 0.5, 0.7], 'LineWidth', 1.2);
h2 = plot3(nan, nan, nan, '*', 'Color', 'k');
h3 = plot3(nan, nan, nan, 'og', 'MarkerFaceColor', 'g');
h4 = plot3(nan, nan, nan, 'om', 'MarkerFaceColor', 'm');
h5 = plot3(nan, nan, nan, 'oy', 'MarkerFaceColor', 'y');

sys = @(t,x) [x(1)*(1-x(1)) - A12*x(1)*x(2) - A13*x(1)*x(3);
              R2*x(2)*(1-x(2)) - A21*x(1)*x(2);
              (R3*x(1)*x(3))/(x(1)+k3) - A31*x(1)*x(3) - D3*x(3)];

% Simulate trajectories from a grid of initial populations.
[x_start, y_start, z_start] = meshgrid([0.05, 0.5], [0.1, 0.9], [0.1, 0.5]);
for k = 1:numel(x_start)
    x0 = [x_start(k); y_start(k); z_start(k)];
    [~, traj] = ode45(sys, [0 350], x0);
    
    plot3(x0(1), x0(2), x0(3), '*', 'Color', 'k', 'MarkerSize', 5);
    plot3(traj(:,1), traj(:,2), traj(:,3), 'Color', [0, 0.25, 0.5, 0.7], 'LineWidth', 1.2);
end

% Collect boundary and coexistence equilibrium candidates.
Eq = [0; 0; 0]; Names = {'P_0'};
Eq = [Eq, [0; 1; 0]]; Names{end+1} = 'P_H';
Eq = [Eq, [1; 0; 0]]; Names{end+1} = 'P_1';

coef_a = A31; 
coef_b = A31*k3 + D3 - R3; 
coef_c = D3*k3;
delta = coef_b^2 - 4*coef_a*coef_c;

if delta >= 0
    roots_x1 = [(-coef_b - sqrt(delta))/(2*coef_a); (-coef_b + sqrt(delta))/(2*coef_a)];
    for x1 = roots_x1'
        if x1 > 0 && x1 < 1.1
            % Check the P13 equilibrium
            x3 = (1 - x1)/A13;
            if x3 >= 0
                Eq = [Eq, [x1; 0; x3]]; Names{end+1} = 'P_{13}';
            end
            % Check the coexistence equilibrium P3
            x2 = 1 - (A21/R2)*x1;
            if x2 > 0
                x3_coex = (1 - x1 - A12*x2)/A13;
                if x3_coex >= 0
                    Eq = [Eq, [x1; x2; x3_coex]]; Names{end+1} = 'P_3';
                end
            end
        end
    end
end

for k = 1:size(Eq, 2)
    p = Eq(:, k);
    
    J = [1 - 2*p(1) - A12*p(2) - A13*p(3),   -A12*p(1),                      -A13*p(1);
         -A21*p(2),                          R2*(1 - 2*p(2)) - A21*p(1),     0;
         (R3*k3*p(3))/((p(1)+k3)^2) - A31*p(3), 0, (R3*p(1))/(p(1)+k3) - A31*p(1) - D3];
     
    % Use the local Jacobian eigenvalues to choose the equilibrium marker.
    vals = real(eig(J));
    
    if any(vals > 1e-5) && any(vals < -1e-5)
        color = 'm';
    elseif any(abs(vals) < 1e-5)
        color = 'y'; 
    else
        color = 'g';
    end
    
    plot3(p(1), p(2), p(3), 'o', 'MarkerSize', 10, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k');
    text(p(1), p(2), p(3)+0.07, Names{k}, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end

legend([h1, h2, h3, h4, h5], {'Trajectories', 'Initial points', 'Stable sink', 'Saddle point', 'Non-hyperbolic'}, ...
    'Location', 'northeast');
