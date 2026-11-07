%% Symbolic Lie derivatives for the tumor-cell output
clear all; clc;

syms x1 x2 x3 real
syms A12 A13 A21 A31 R2 R3 k3 D3 real

x = [x1; x2; x3];

% 1. Drift vector field f(x)
f = [
    x1*(1 - x1) - A12*x1*x2 - A13*x1*x3;        % f1
    R2*x2*(1 - x2) - A21*x1*x2;                 % f2
    (R3*x1*x3)/(x1 + k3) - A31*x1*x3 - D3*x3    % f3
];

% 2. Lie derivatives for output y = x1
% h(x) = x1
Lf_h = f(1); % First derivative (x1_dot)

% L_f^2 h: second derivative along f
Lf2_h = jacobian(Lf_h, x) * f;
Lf2_h = simplify(Lf2_h);

% 3. Display the expression
fprintf('--- SYMBOLIC EXPRESSION FOR THE SIMULINK BLOCK ---\n\n');
fprintf('Lf2_h = %s;\n', string(Lf2_h));
