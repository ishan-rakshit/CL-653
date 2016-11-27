function [params,mats,init,x,y] = true_states(params,mats,init)

% This case study is taken from Farina et al. (2002) and more details can
% be found in this reference.
%     s = rng; % the new random number generator seed function
    randn('state',0) % initializing the random number generator
    
    % PROCESS MODEL
    
    % State vector (x1 x3 - positions in x and y, x2 x4 - velocities in x
    % and y)
    x = zeros(params.dim,params.N); % true states captured in 4x1 vectors arranged for all N measurements
    y = zeros(2,params.N);
    
    % Initial conditions for true states generation
    x(:,1) = [init.x0 -2.2552e+03 init.y0 -397.6543]';
%     [~,~,rho] = return_rho(x(3,1));
    v = mvnrnd(zeros(1,2),mats.R)';
    r = sqrt(init.y0^2 + init.x0^2) + v(1,1);
    eps = init.gamma0;
    sigma_d = mats.R(1,1)^2*(cos(eps))^2 + r^2*mats.R(2,2)*(sin(eps))^2;
    sigma_h = mats.R(1,1)^2*(sin(eps))^2 + r^2*mats.R(2,2)*(cos(eps))^2;
    sigma_dh = (mats.R(1,1)^2 - r^2*mats.R(2,2)^2)*sin(eps)*cos(eps);
    sigma_d = sqrt(sigma_d);
    sigma_h = sqrt(sigma_h);
    sigma_dh = sqrt(sigma_dh);
    init.P = [sigma_d -sigma_d/params.Ts sigma_dh -sigma_dh/params.Ts; -sigma_d/params.Ts 2*sigma_d/params.Ts^2 -sigma_dh/params.Ts 2*sigma_dh/params.Ts^2; sigma_dh -sigma_dh/params.Ts sigma_h -sigma_h/params.Ts; -sigma_dh/params.Ts 2*sigma_dh/params.Ts^2 -sigma_h/params.Ts 2*sigma_h/params.Ts^2]; %covariance matrix

    % Process noise
    q = 1; % parameter related to process noise intensity
    mats.theta_mat = q*[params.Ts^3/3 params.Ts^2/2; params.Ts^2/2 params.Ts];
    mats.Q_d = [mats.theta_mat zeros(size(mats.theta_mat)); zeros(size(mats.theta_mat)) mats.theta_mat]; % covariance matrix of the process noise
    w = mvnrnd(zeros(1,4),mats.Q_d,params.N)';
    
    % generating true states from the initial condition given
    for k = 1:params.N-1
        f = return_f(x(:,k),params);
        x(:,k+1) = mats.Phi*x(:,k) + mats.G*f + mats.G*[0; -params.g] + w(:,k);
        y(:,k) = mats.C*x(:,k) + v; 
        r = sqrt(y(1,k)^2 + y(2,k)^2); 
        eps = atan2(y(2,k),y(1,k));
        sigma_d = mats.R(1,1)^2*(cos(eps))^2 + r^2*mats.R(2,2)*(sin(eps))^2;
        sigma_h = mats.R(1,1)^2*(sin(eps))^2 + r^2*mats.R(2,2)*(cos(eps))^2;
        sigma_dh = (mats.R(1,1)^2 - r^2*mats.R(2,2)^2)*sin(eps)*cos(eps);
        sigma_d = sqrt(sigma_d);
        sigma_h = sqrt(sigma_h);
        sigma_dh = sqrt(sigma_dh);
        Rk = [sigma_d sigma_dh; sigma_dh sigma_h];
        v = mvnrnd(zeros(1,2),Rk)';
    end
    v = mvnrnd(zeros(1,2),Rk)';
    y(:,params.N) = mats.C*x(:,params.N)+v;
    
    % Let the plotting begin!
%     figure
%     subplot(2,1,1)
%     plot(x(1,:),x(3,:),'k')
%     title('position')
%     subplot(2,1,2)
%     plot((x(2,:).^2+x(4,:).^2),'b')
%     title('velocity')
end