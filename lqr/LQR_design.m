%% LQR design near the tumor-free equilibrium
clc; clear;

p.R2 = 0.6; 
p.D3 = 0.5; 
p.A12 = 1.0; 
p.A13 = 2.5; 
p.A21 = 1.5; 
p.A31 = 0.2; 
p.k3 = 1.0;
p.R3 = 4.5; 

% Select an operating point close to the tumor-free equilibrium.
epsilon = 0.01;
x_lin = [epsilon; 1; 0];
x1=x_lin(1); x2=x_lin(2); x3=x_lin(3);

% Linearized state matrix and control-input direction.
A = [(1 - 2*x1 - p.A12*x2 - p.A13*x3),   -p.A12*x1,                     -p.A13*x1;
     -p.A21*x2,                          p.R2*(1 - 2*x2) - p.A21*x1,     0;
     (p.R3*x3*p.k3)/(x1 + p.k3)^2 - p.A31*x3, 0,                        (p.R3*x1)/(x1 + p.k3) - p.A31*x1 - p.D3];

B = [0; 0; 1];

if rank(ctrb(A,B)) == 3
    disp('The system is controllable at P_epsilon.');
else
    error('The system is not controllable at P_epsilon.');
end

% State and input weights in the quadratic cost.
Q = diag([100, 1, 1]);
R = 1;

K = lqr(A, B, Q, R);

fprintf('\nLQR feedback gain K:\n');
disp(K);

% Desired population vector for the closed-loop simulation.
x_target = [0; 1; 0];
