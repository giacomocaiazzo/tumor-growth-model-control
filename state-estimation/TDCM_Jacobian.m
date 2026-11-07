function Ad = Jacobian(x, u)
% State-transition Jacobian used by the Simulink extended Kalman filter.
   
    A12 = 0.1; A13 = 0.1;
    R2 = 1.5;  A21 = 0.1;
    R3 = 4.5; k3 = 1.0; A31 = 0.2; D3 = 0.5;
    
    dt = 0.001;
    x1 = x(1); x2 = x(2); x3 = x(3);

    A11 = 1 - 2*x1 - A12*x2 - A13*x3; 
    A12 = -A12 * x1;                 
    A13 = -A13 * x1;                 
    A21 = -A21 * x2;                  
    A22 = R2*(1 - 2*x2) - A21*x1;    
    A23 = 0;                          
    d_rec_dx1 = (R3*x3 * (x1+k3) - R3*x1*x3 * 1) / (x1 + k3)^2; 
    A31 = (R3*x3*k3)/(x1 + k3)^2 - A31*x3;
    A32 = 0;
    A33 = (R3*x1)/(x1 + k3) - A31*x1 - D3;
    
    A = [A11, A12, A13; 
         A21, A22, A23; 
         A31, A32, A33];

    % First-order discrete-time approximation for the EKF prediction.
    Ad = eye(3) + A * dt;
end
