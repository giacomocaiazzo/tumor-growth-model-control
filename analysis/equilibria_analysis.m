%% Open-loop equilibrium and stability analysis
% Solve f(x) = 0 and evaluate the Jacobian eigenvalues.

% Run tumor_model_parameters.m first to define the parameter structure p.
if ~exist('p', 'var')
    error('Run tumor_model_parameters.m first.');
end

% Autonomous system (u = 0) for the equilibrium solver.
% x = [x1; x2; x3]
F = @(x) [
    x(1)*(1 - x(1)) - p.A12*x(1)*x(2) - p.A13*x(1)*x(3);        % dot_x1
    p.R2*x(2)*(1 - x(2)) - p.A21*x(1)*x(2);                     % dot_x2
    (p.R3*x(1)*x(3))/(x(1) + p.k3) - p.A31*x(1)*x(3) - p.D3*x(3) % dot_x3
];

% --- 1. Equilibrium search ---
options = optimoptions('fsolve', 'Display', 'off');

% Use multiple initial guesses to search for distinct equilibria.
x0_guesses = [0 0 0; 0 1 0; 1 0 0; 0.5 0.5 0.5; 1 1 1; 0.1 0.6 0.1];
equilibria = [];

for i = 1:size(x0_guesses, 1)
    [x_eq, fval, exitflag] = fsolve(F, x0_guesses(i,:)', options);
    if exitflag > 0
        % Keep nonnegative solutions and discard duplicates.
        if all(x_eq >= -1e-6) % Tolerance for numerical zeros
            x_eq = round(x_eq, 4); % Round to identify duplicate solutions
            if isempty(equilibria) || ~ismember(x_eq', equilibria, 'rows')
                equilibria = [equilibria; x_eq'];
            end
        end
    end
end

% --- 2. Jacobian-based stability analysis ---
disp('--------------------------------------------------');
disp('OPEN-LOOP EQUILIBRIUM ANALYSIS');
disp('--------------------------------------------------');

for i = 1:size(equilibria, 1)
    eq = equilibria(i, :)';
    
    % Evaluate the analytical Jacobian at the equilibrium.
    % J_11 = df1/dx1 ...
    x1 = eq(1); x2 = eq(2); x3 = eq(3);
    
    % Partial derivatives
    df1_dx1 = 1 - 2*x1 - p.A12*x2 - p.A13*x3;
    df1_dx2 = -p.A12*x1;
    df1_dx3 = -p.A13*x1;
    
    df2_dx1 = -p.A21*x2;
    df2_dx2 = p.R2 - 2*p.R2*x2 - p.A21*x1;
    df2_dx3 = 0;
    
    % Recruitment derivative using the quotient rule
    % N = R3*x1*x3, D = x1+k3
    % d/dx1 = (R3*x3 * (x1+k3) - R3*x1*x3 * 1) / (x1+k3)^2 
    %       = (R3*x3*k3) / (x1+k3)^2
    dRec_dx1 = (p.R3 * x3 * p.k3) / (x1 + p.k3)^2;
    df3_dx1 = dRec_dx1 - p.A31*x3;
    
    df3_dx2 = 0;
    
    % d/dx3 of (R3*x1*x3)/(x1+k3) is (R3*x1)/(x1+k3)
    dRec_dx3 = (p.R3 * x1) / (x1 + p.k3);
    df3_dx3 = dRec_dx3 - p.A31*x1 - p.D3;
    
    J = [df1_dx1, df1_dx2, df1_dx3;
         df2_dx1, df2_dx2, df2_dx3;
         df3_dx1, df3_dx2, df3_dx3];
     
    eigs = eig(J);
    
    % Display the results
    fprintf('Equilibrium E%d: [%.4f, %.4f, %.4f]\n', i-1, x1, x2, x3);
    disp('Eigenvalues:');
    disp(eigs);
    
    if all(real(eigs) < 0)
        disp('--> STABLE (sink)');
    elseif all(real(eigs) > 0)
        disp('--> UNSTABLE (source)');
    else
        disp('--> UNSTABLE (saddle)');
    end
    disp(' ');
end
