% TRACKING A BALLISTIC TARGET
% For KF implementation --
% find a steady state around which to obtain a linearized single model. 
% The steady state can be either obtained analytically or by noise-free 
% simulation of the given model

function [x_hat_updated,error_KF,rho_predicted,rho_updated] = filter_kalman(init,params,mats,x,y)
% We have to linearize the system about such a point that is the stationary 
% state for a general discrete system xn+1 = f(xn) is a x1, such that x1 = f(x1)
    x_steady_state = x(:,1);
    x_difference = zeros(4,params.N);
    y_steady_state = mats.C*x_steady_state; 
    F = return_jacobian(params,x(:,1)); % linearization about the same point
    x_hat_updated = zeros(4,params.N);
    x_hat_updated(:,1) = x(:,1);
    x_difference(:,1) = x_hat_updated(:,1) - x_steady_state;
    error_KF = zeros(2,params.N-1);
    rho_predicted = zeros(1,params.N);
    rho_updated = zeros(1,params.N);
    P_hat_updated = init.P;
    rho_predicted(1) = max(eig(P_hat_updated));
    rho_updated(1) = max(eig(P_hat_updated));
    
    for k = 1:params.N-1
        r = sqrt(y(2,k)^2 + y(1,k)^2); 
        eps_KF = atan2(y(2,k),y(1,k));
        sigma_d = mats.R(1,1)^2*(cos(eps_KF))^2 + r^2*mats.R(2,2)*(sin(eps_KF))^2;
        sigma_h = mats.R(1,1)^2*(sin(eps_KF))^2 + r^2*mats.R(2,2)*(cos(eps_KF))^2;
        sigma_dh = (mats.R(1,1)^2 - r^2*mats.R(2,2)^2)*sin(eps_KF)*cos(eps_KF);
        sigma_d = sqrt(sigma_d);
        sigma_h = sqrt(sigma_h);
        sigma_dh = sqrt(sigma_dh);
        Rk = [sigma_d sigma_dh; sigma_dh sigma_h];
        x_hat_predicted = (mats.Phi + mats.G*F)*x_difference(:,k)+ mats.G*[0; -params.g];
        P_hat_predicted = (mats.Phi + mats.G*F)*P_hat_updated*(mats.Phi + mats.G*F)' + mats.Q_d;
        rho_predicted(k+1) = max(eig(P_hat_predicted));
        L = P_hat_predicted*mats.C'/(mats.C*P_hat_predicted*mats.C' + Rk);
        error_KF(:,k) = y(:,k)- mats.C*x_hat_predicted - y_steady_state;
        x_difference(:,k+1) = x_hat_predicted + L*error_KF(:,k);
        x_hat_updated(:,k+1) = x_steady_state + x_difference(:,k+1);
        P_hat_updated = (eye(4)-L*mats.C)*P_hat_predicted;
        rho_updated(k+1) = max(eig(P_hat_updated));
    end
end