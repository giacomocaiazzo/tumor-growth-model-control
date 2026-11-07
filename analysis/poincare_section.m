%% Poincare section of the tumor-growth model
clear; clc; close all;

% Parameter set for the chaotic regime
p.R2=0.6; p.R3=4.5; p.D3=0.5; p.A12=1.0; p.A13=2.5; p.A21=1.5; p.A31=0.2; p.k3=1.0; u=0;

% Use a long trajectory to collect section crossings.
T_max = 5000; 
x0 = [0.1; 0.1; 0.1];

% Dimensionless model equations
cancer_ode = @(t, x) [
    x(1)*(1 - x(1)) - p.A12*x(1)*x(2) - p.A13*x(1)*x(3);
    p.R2*x(2)*(1 - x(2)) - p.A21*x(1)*x(2);
    (p.R3*x(1)*x(3))/(x(1) + p.k3) - p.A31*x(1)*x(3) - p.D3*x(3) + u
];

disp('Simulating the open-loop trajectory...');
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
[t, x] = ode45(cancer_ode, [0 T_max], x0, options);

% --- Poincare section ---
% Section plane: x2 = 0.6
% Keep crossings in the decreasing x2 direction.
plane_val = 0.6;
poincare_points = [];

% Discard the first 30% of the trajectory as a transient.
start_idx = round(length(x) * 0.3);

for i = start_idx : length(x)-1
    % Detect crossings from above to below the section plane.
    if (x(i, 2) > plane_val) && (x(i+1, 2) < plane_val)
        % Estimate the crossing coordinates by linear interpolation.
        delta = (plane_val - x(i, 2)) / (x(i+1, 2) - x(i, 2));
        
        x1_cross = x(i, 1) + delta * (x(i+1, 1) - x(i, 1));
        x3_cross = x(i, 3) + delta * (x(i+1, 3) - x(i, 3));
        
        poincare_points = [poincare_points; x1_cross, x3_cross];
    end
end

% --- Plot ---
figure('Color','w');
plot(poincare_points(:,1), poincare_points(:,2), 'k.', 'MarkerSize', 4);
xlabel('x_1 (Tumor Cells)');
ylabel('x_3 (Resting Cells)');
title(['Poincare section at x_2 = ', num2str(plane_val)]);
grid on;

% Display a qualitative classification of the sampled section.
if size(poincare_points,1) < 10
    disp('Qualitative classification: possible limit cycle (few sampled points)');
else
    disp('Qualitative classification: possible torus or strange attractor (inspect the section)');
end
