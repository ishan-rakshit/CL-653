% TRACKING A BALLISTIC TARGET
% For KF implementation --
% find a steady state around which to obtain a linearized single model. 
% The steady state can be either obtained analytically or by noise-free 
% simulation of the given model

function filter_kalman(params,mats,x)
    x_hat_predicted = zeros(4,params.N);
    x_hat_updated = zeros(4,params.N);
    P_hat_predicted = zeros(4,4);
    P_hat_updated = zeros(4,4);
    
    % Initial conditions for the filter are specified
    x1 = x(1,1); x2 = x(2,1); x3 = x(3,1); x4 = x(4,1);
%     x_hat_updated(:,1) = [x2 (x2-x1)/params.Ts y2 (y2-y1)/params.Ts]';
    x_hat_updated(:,1) = [x(1,2) (x(1,2)-x1)/params.Ts x(3,2) (x(3,2)-x3)/params.Ts]';
    % Initial covariance matrix
    sigma_r = 100;
    sigma_epsilon = 0.017;
    r = sqrt(x1^2 + x3^2);
    epsilon = atan(x3/x1);
    sigma_d = sqrt(sigma_r^2 * cos(epsilon)^2 + r^2 * sigma_epsilon^2 * sin(epsilon)^2);
    sigma_h = sqrt(sigma_r^2 * sin(epsilon)^2 + r^2 * sigma_epsilon^2 * cos(epsilon)^2);
    sigma_dh = (sigma_r^2  - r^2 * sigma_epsilon^2) * sin(epsilon)*cos(epsilon);
    P_hat_updated = [sigma_d^2  -sigma_d^2/params.Ts    sigma_dh    -sigma_dh/params.Ts;
        -sigma_d^2/params.Ts    2*sigma_d^2/(params.Ts^2)   -sigma_dh/params.Ts   2*sigma_dh/(params.Ts^2);
        sigma_dh    -sigma_dh/params.Ts sigma_h^2   -sigma_h^2/params.Ts;
        -sigma_dh/params.Ts     2*sigma_dh/params.Ts^2  -sigma_h^2/params.Ts    2*sigma_h^2/params.Ts^2];
    
    for k = 2:params.N
        % prediction step
        x_hat_predicted(:,k-1) = mats.Phi*x_hat_updated(:,k-1);
        P_hat_predicted = mats.Phi*P_hat_updated*mats.Phi';
        % kalman gain computation
        L_star = P_hat_predicted*C'/(C*P_hat_predicted*C' + R);
        % update step
        e = y(k) - C*x_hat_predicted;
        P_hat_updated = (eye(size(L_star*C)) - L_star*C)*P_hat_predicted;
        x_hat_updated(:,k) = x_hat_predicted(:,k-1) + L_star*e;
    end
end