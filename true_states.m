function [params,mats,x] = true_states

% This case study is taken from Farina et al. (2002) and more details can
% be found in this reference.
    s = rng;
%     randn('state',5) % initializing the random number generator
    
    % PROCESS MODEL
    
    % Given constant parameters for the model
    params.beta = 40000; % kg/(ms^2) - ballistic coefficient
    params.g = 9.81; % m/s^2 - acceleration due to gravity
    params.Ts = 2; % s - sampling interval
    params.q = 1; % m^2/s^3 - process noise intensity
    params.N = 60; % time steps being considered (NOT SPECIFIED IN PROBLEM)
    params.dim = 4; % dimension of the state vector
    
    % State vector (x1 x3 - positions in x and y, x2 x4 - velocities in x
    % and y)
    x = zeros(params.dim,params.N); % true states captured in 4x1 vectors arranged for all N measurements
    
    % Initial conditions for true states generation
    x(:,1) = [232000 2.2552e+03 88000 397.6543]';
    
    % Process noise
    q = 1; % parameter related to process noise intensity
    mats.theta_mat = q*[params.Ts^3/3 params.Ts^2/2; params.Ts^2/2 params.Ts];
    mats.Q_d = [mats.theta_mat zeros(size(mats.theta_mat)); zeros(size(mats.theta_mat)) mats.theta_mat]; % covariance matrix of the process noise
    w = mvnrnd(zeros(1,4),mats.Q_d,params.N)';
    
    % defining matrices
    mats.Phi = [1 params.Ts 0 0; 0 1 0 0; 0 0 1 params.Ts; 0 0 0 1];
    mats.G = [params.Ts^2/2 0; params.Ts 0; 0 params.Ts^2/2; 0 params.Ts];
    mats.C = [1 0 0 0; 0 0 1 0];
    % generating true states from the initial condition given
    for k = 1:params.N-1
        x(:,k+1) = return_psi(x(:,k),mats,params) + mats.G*[0 -params.g]' + w(:,k);
    end
    
    rng(s);
end