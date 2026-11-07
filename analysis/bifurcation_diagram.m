%% Period-doubling bifurcation diagram
clear; clc; close all;

% Model parameters
p.R2=0.6; p.D3=0.5; p.A12=1.0; p.A13=2.5; p.A21=1.5; p.A31=0.2; p.k3=1.0;

% Sweep values, sampled peaks and initial state
R3_vec = 2.5:0.002:4.5;
R3_plot = []; Peaks_plot = [];
x0 = [0.2; 0.6; 0.2];

for i = 1:length(R3_vec)
    p.R3 = R3_vec(i);
    [~, x_tr] = ode45(@(t,x) sys_ode(t,x,p), [0 1500], x0); % Discard a long transient before collecting peaks
    x0_reg = x_tr(end, :)';
    [~, x_reg] = ode45(@(t,x) sys_ode(t,x,p), [0 400], x0_reg);
    tumor = x_reg(:,1);
    [pks, ~] = findpeaks(tumor);
    if length(pks) > 40, pks = pks(end-35:end); end
    R3_plot = [R3_plot; ones(length(pks),1) * p.R3];
    Peaks_plot = [Peaks_plot; pks];
    % Continue the sweep from the final transient state.
    x0 = x0_reg;
end

% --- Bifurcation plot ---
figure('Color', 'w', 'Name', 'Period-Doubling Bifurcation Diagram');
scatter(R3_plot, Peaks_plot, 1, 'k', 'filled', 'MarkerEdgeAlpha', 0.1);
hold on;

% Reference bifurcation values used in the report
R_hopf = 2.686;
R_flip1 = 3.293; % First period doubling
R_flip2 = 3.447; % Second period doubling

xline(R_hopf, 'r--', 'Hopf Bif.', 'LabelVerticalAlignment', 'bottom');
xline(R_flip1, 'b--', '1st period doubling (T -> 2T)', 'LabelVerticalAlignment', 'bottom');
xline(R_flip2, 'm--', '2nd period doubling (2T -> 4T)', 'LabelVerticalAlignment', 'bottom');

title('Bifurcation Diagram: Period-Doubling Cascade');
xlabel('Recruitment Rate R_3'); ylabel('Tumor Peaks x_1');
grid on; axis([2.5 4.5 0 1]);

% Uncontrolled dynamics with states [tumor; hunting cells; resting cells].
function dx = sys_ode(~, x, p)
    dx = [x(1)*(1-x(1)) - p.A12*x(1)*x(2) - p.A13*x(1)*x(3);
          p.R2*x(2)*(1-x(2)) - p.A21*x(1)*x(2);
          (p.R3*x(1)*x(3))/(x(1)+p.k3) - p.A31*x(1)*x(3) - p.D3*x(3)];
end
