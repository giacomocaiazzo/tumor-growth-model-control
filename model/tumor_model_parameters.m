%% Tumor-growth model parameters and initial conditions
% State order: tumor, hunting cells, resting cells.
ip=0;
p = struct();

p.R2  = 0.6;    
p.R3  = 2.69;    
p.D3  = 0.5;   
p.A12 = 1.0;    
p.A13 = 2.5;    
p.A21 = 1.5;    
p.A31 = 0.2;    
p.k3  = 1.0;    

% Create the parameter bus used by the Simulink models.
Simulink.Bus.createObject(p);

% Initial populations in the same order as the state vector.
x1_0 = 0.1; 
x2_0 = 0.1; 
x3_0 = 0.1; 
x0 = [x1_0; x2_0; x3_0];


