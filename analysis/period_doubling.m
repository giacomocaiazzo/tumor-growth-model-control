%% Period-doubling detection
clear; clc; close all;

% Model parameters
p.R2=0.6; p.D3=0.5; p.A12=1.0; p.A13=2.5; p.A21=1.5; p.A31=0.2; p.k3=1.0; u=0;

% --- Parameter sweep ---
% Start beyond the Hopf bifurcation near R3 = 2.68.
r3_start = 2.8; 
r3_end = 4.0;
r3_vec = r3_start:0.002:r3_end;

x0 = [0.2; 0.6; 0.2]; % Initial condition
current_period = 1;   % Initial period-1 limit-cycle assumption
bifurcations = [];

fprintf('Scanning period doublings (range %.2f - %.2f)...\n', r3_start, r3_end);

for i = 1:length(r3_vec)
    p.R3 = r3_vec(i);
    
    % 1. Discard a long transient before sampling the response.
    [~, x_tr] = ode45(@(t,x) sys_ode(t,x,p), [0 1000], x0);
    % Initialize the sampling interval from the final transient state.
    x0 = x_tr(end, :)'; 
    
    % 2. Sample the steady-state response.
    [t, x] = ode45(@(t,x) sys_ode(t,x,p), [0 300], x0);
    signal = x(:,1); % Tumor population
    
    % 3. Peak analysis
    [pks, ~] = findpeaks(signal);
    
    % Screen out stationary responses.
    % If the oscillation decays between the first and last peak,
    % skip the response as a fixed point or a remaining transient.
    if isempty(pks) || (max(pks) - min(pks) < 1e-4)
        % Skip the fixed-point response.
        continue; 
    end
    
    % Analyze the remaining oscillatory response.
    % Count distinct peak levels within a tolerance.
    % Use the tolerance below to distinguish period doublings.
    unique_pks = uniquetol(pks, 1e-2); 
    n_period = length(unique_pks);
    
    % 4. Detect period doublings: 1 -> 2, 2 -> 4, etc.
    if n_period == current_period * 2
        fprintf('>>> FLIP BIFURCATION! R3 = %.4f | Period %d -> %d\n', ...
            p.R3, current_period, n_period);
        bifurcations = [bifurcations; p.R3, current_period, n_period];
        current_period = n_period;
        
        % Stop after detecting a period of at least eight.
        if current_period >= 8
            fprintf('    (Period-doubling cascade detected)\n');
            break;
        end
    elseif n_period > current_period && n_period ~= current_period * 2
        % Other periods, such as three or five, may indicate periodic windows.
        % Do not record these as standard period doublings.
    end
end

% --- Model equations ---
function dx = sys_ode(~, x, p)
    x1=x(1); x2=x(2); x3=x(3);
    dx = zeros(3,1);
    dx(1) = x1*(1-x1) - p.A12*x1*x2 - p.A13*x1*x3;
    dx(2) = p.R2*x2*(1-x2) - p.A21*x1*x2;
    dx(3) = (p.R3*x1*x3)/(x1+p.k3) - p.A31*x1*x3 - p.D3*x3;
end
