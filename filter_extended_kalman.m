% TRACKING A BALLISTIC TARGET
% For EKF implementation --

function [x_hat_updated,error_EKF,rho_predicted,rho_updated] = filter_extended_kalman(init,params,mats,x,y)
    x_hat_updated = zeros(4,params.N);
    x_hat_updated(:,1) = x(:,1);
    error_EKF = zeros(2,params.N-1);
    rho_predicted = zeros(1,params.N);
    rho_updated = zeros(1,params.N);
    P_hat_updated = init.P;
    rho_predicted(1) = max(eig(P_hat_updated));
    rho_updated(1) = max(eig(P_hat_updated));

    for k = 1:params.N-1
        F = return_jacobian(params,x_hat_updated(:,k));
        f = return_f(x_hat_updated(:,k),params);
        x_hat_predicted = mats.Phi*x_hat_updated(:,k) + mats.G*f + mats.G*[0; -params.g] ;
    %     x_pred = (phi + G*F)*x_EKF(:,i) + G*[0;-g];
        P_hat_predicted = (mats.Phi + mats.G*F)*P_hat_updated*(mats.Phi + mats.G*F)' + mats.Q_d;
        rho_predicted(k+1) = max(eig(P_hat_predicted));
        r = sqrt(y(2,k)^2 + y(1,k)^2); 
        eps_EKF = atan2(y(2,k),y(1,k));
        
        sigma_d = mats.R(1,1)^2*(cos(eps_EKF))^2 + r^2*mats.R(2,2)*(sin(eps_EKF))^2;
        sigma_h = mats.R(1,1)^2*(sin(eps_EKF))^2 + r^2*mats.R(2,2)*(cos(eps_EKF))^2;
        sigma_dh = (mats.R(1,1)^2 - r^2*mats.R(2,2)^2)*sin(eps_EKF)*cos(eps_EKF);
        sigma_d = sqrt(sigma_d);
        sigma_h = sqrt(sigma_h);
        sigma_dh = sqrt(sigma_dh);
        Rk = [sigma_d sigma_dh; sigma_dh sigma_h];
        L = P_hat_predicted*mats.C'/(mats.C*P_hat_predicted*mats.C' + Rk);
        error_EKF(:,k) = y(:,k) - mats.C*x_hat_predicted; 
        x_hat_updated(:,k+1) = x_hat_predicted + L*error_EKF(:,k);
        P_hat_updated = (eye(4) - L*mats.C)*P_hat_updated;
        rho_updated(k+1) = max(eig(P_hat_updated));
    end




%     x_hat_predicted = zeros(params.dim,params.N);
%     x_hat_updated = zeros(params.dim,params.N);
% %     P_hat_predicted = zeros(params.dim,params.dim);
% %     P_hat_updated = zeros(params.dim,params.dim);
%     
%     % Initial conditions for the filter are specified
%     x1 = x(1,1); 
% %     x2 = x(2,1); 
%     x3 = x(3,1); 
% %     x4 = x(4,1);
%     
% %     x_hat_updated(:,1) = [x2 (x2-x1)/params.Ts y2 (y2-y1)/params.Ts]';
%     x_hat_updated(:,1) = [x(1,2) (x(1,2)-x1)/params.Ts x(3,2) (x(3,2)-x3)/params.Ts]';
%     % Initial covariance matrix
%     sigma_r = 100;
%     sigma_epsilon = 0.017;
%     r = sqrt(x1^2 + x3^2);
%     epsilon = atan(x3/x1);
%     sigma_d = sqrt(sigma_r^2 * cos(epsilon)^2 + r^2 * sigma_epsilon^2 * sin(epsilon)^2);
%     sigma_h = sqrt(sigma_r^2 * sin(epsilon)^2 + r^2 * sigma_epsilon^2 * cos(epsilon)^2);
%     sigma_dh = (sigma_r^2  - r^2 * sigma_epsilon^2) * sin(epsilon)*cos(epsilon);
%     P_hat_updated = [sigma_d^2  -sigma_d^2/params.Ts    sigma_dh    -sigma_dh/params.Ts;
%         -sigma_d^2/params.Ts    2*sigma_d^2/(params.Ts^2)   -sigma_dh/params.Ts   2*sigma_dh/(params.Ts^2);
%         sigma_dh    -sigma_dh/params.Ts sigma_h^2   -sigma_h^2/params.Ts;
%         -sigma_dh/params.Ts     2*sigma_dh/params.Ts^2  -sigma_h^2/params.Ts    2*sigma_h^2/params.Ts^2];
%     
%     for k = 2:params.N
%         % compute jacobian for this step
%         Psi_k = return_jacobian(params, mats, x_hat_updated(:,k-1));
%         % prediction step
%         x_hat_predicted(:,k-1) = Psi_k*x_hat_updated(:,k-1) + mats.G * [0 -params.g]';
%         P_hat_predicted = Psi_k*P_hat_updated*Psi_k' + mats.Q_d;
%         % kalman gain computation
%         L_star = P_hat_predicted*mats.C'/(mats.C*P_hat_predicted*mats.C' + mats.R);
%         % update step
%         e = y(:,k) - mats.C*x_hat_predicted(:,k-1);
%         P_hat_updated = (eye(size(L_star*mats.C)) - L_star*mats.C)*P_hat_predicted;
%         x_hat_updated(:,k) = x_hat_predicted(:,k-1) + L_star*e;
%     end
end