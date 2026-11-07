%% Symbolic Hopf bifurcation analysis
clear; clc;

% 1. Symbolic variables
syms x1 x2 x3 R3 real
% Fixed model parameters
R2=0.6; D3=0.5; A12=1.0; A13=2.5; A21=1.5; A31=0.2; k3=1.0;

% 2. Symbolic Jacobian at the coexistence equilibrium
% Recruitment coupling derivative
d_recruitment = (R3*k3*x3)/((x1 + k3)^2) - A31*x3; % Entry J(3,1)

J = [ -x1,      -A12*x1,          -A13*x1;
      -A21*x2,  -R2*x2,           0;
      d_recruitment, 0,           0 ]; 
% The equilibrium relations give J(3,3) = 0 and J(1,1) = -x1.

% 3. Characteristic polynomial coefficients
% P(lambda) = lambda^3 + a1*lambda^2 + a2*lambda + a3
p_char = charpoly(J);
a1 = p_char(2);
a2 = p_char(3);
a3 = p_char(4);

% 4. Critical Routh-Hurwitz condition
Hopf_Condition = a1 * a2 - a3;

fprintf('Symbolic Hopf condition (a1*a2 - a3 = 0):\n');
pretty(simplify(Hopf_Condition));

% 5. Substitute the equilibrium relations
% Express x2 and x3 in terms of x1 and R3.
x2_eq = 1 - (A21/R2)*x1;
x3_eq = (1 - x1 - A12*x2_eq)/A13;

% Substitute x2 and x3 into the Hopf condition.
Hopf_x1_R3 = subs(Hopf_Condition, {x2, x3}, {x2_eq, x3_eq});

% Solve two equations for x1 and R3:
% 1) Hopf condition: Hopf_x1_R3 = 0
% 2) Equilibrium condition: A*x1^2 + B*x1 + C = 0
A_quad = A31;
B_quad = A31*k3 + D3 - R3;
C_quad = D3*k3;
Eq_Equilibrio = A_quad*x1^2 + B_quad*x1 + C_quad;

fprintf('\nSolving for the critical R3 value...\n');

% Solve the nonlinear algebraic system numerically.
% Search for positive x1 and R3.
sol = vpasolve([Hopf_x1_R3 == 0, Eq_Equilibrio == 0], [x1, R3], [0 1; 2 3]);

x1_hopf = double(sol.x1);
R3_hopf = double(sol.R3);

fprintf('\n>>> HOPF ANALYSIS RESULT <<<\n');
fprintf('Hopf bifurcation at R3 = %.6f\n', R3_hopf);
fprintf('Critical tumor coordinate x1 = %.6f\n', x1_hopf);
