function x_next = state_transition(x, u)
% Discrete state transition used by the Simulink extended Kalman filter.

    A12 = 0.1; A13 = 0.1;
    R2 = 1.5;  A21 = 0.1; 
    R3 = 4.5; k3 = 1.0; A31 = 0.2; D3 = 0.5;
    
    dt = 0.001;

    x1 = x(1); x2 = x(2); x3 = x(3);

    dx1 = x1*(1 - x1) - A12*x1*x2 - A13*x1*x3;
    dx2 = R2*x2*(1 - x2) - A21*x1*x2;

 
    % The control input acts on the resting-cell population.
    recruitment = (R3 * x1 * x3) / (x1 + k3);
    dx3 = recruitment - A31*x1*x3 - D3*x3 + u;

    % Forward Euler prediction over one estimator sampling interval.
    x_next = x + [dx1; dx2; dx3] * dt;
end
