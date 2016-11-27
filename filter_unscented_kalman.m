% TRACKING A BALLISTIC TARGET
% For UKF implementation --

function [x_hat_updated,error_UKF,rho_predicted,rho_updated] = filter_unscented_kalman(init,params,mats,x,y)
    P_hat_updated = init.P;
    rho_predicted = zeros(1,params.N);
    rho_updated = zeros(1,params.N);
    rho_predicted(1) = max(eig(P_hat_updated));
    rho_updated(1) = max(eig(P_hat_updated));
    error_UKF = zeros(2,params.N-1);
    W = ones(1,10).*1/(10); % to initialize W
    size_x = size(x);
    W(1) = params.kappa/(size_x(1) + params.kappa);
    x_hat_updated = zeros(4,params.N);
    new_P = sqrtm(P_hat_updated.*(size_x(1) + params.kappa));
    x_hat_updated(:,1) = x(:,1); % initialization using the starting values
    zeta = zeros(2,2*size_x(1) + 1);
    for k = 1:params.N-1
        sigma_updated = zeros(4,2*size_x(1)+1);
        sigma_updated(:,1) = x_hat_updated(:,k);
        for j = 2:size_x(1)+1
            sigma_updated(:,j) = x_hat_updated(:,k) + new_P(:,(j-1));
        end
        for j = size_x(1)+2 : 2*size_x(1)+1
            sigma_updated(:,j) = x_hat_updated(:,k) - new_P(:,(j-size_x(1)-1));
        end
        sigma_predicted = zeros(4,2*size_x(1)+1);
        x_hat_predicted = zeros(4,1);
        P_hat_updated = P_hat_updated + mats.Q_d;
        y_predicted = zeros(2,1);
        for j = 1:2*size_x(1) +1
            f = return_f(sigma_updated(:,j),params);
            sigma_predicted(:,j) = mats.Phi*sigma_updated(:,j) + mats.G*f + mats.G*[0; -params.g];
            x_hat_predicted = x_hat_predicted + sigma_predicted(:,j).*W(j);
            P_hat_updated = P_hat_updated + (sigma_predicted(:,j) - x_hat_predicted)*(sigma_predicted(:,j) - x_hat_predicted)'.*W(j);
            zeta(:,j) = mats.C*sigma_predicted(:,j);
            y_predicted = y_predicted + zeta(:,j).*W(j);
        end
        rho_predicted(k+1) = max(eig(P_hat_updated));
        r = sqrt(y(2,k)^2 + y(1,k)^2); 
        eps = atan2(y(2,k),y(1,k));
        sigma_d = mats.R(1,1)^2*(cos(eps))^2 + r^2*mats.R(2,2)*(sin(eps))^2;
        sigma_h = mats.R(1,1)^2*(sin(eps))^2 + r^2*mats.R(2,2)*(cos(eps))^2;
        sigma_dh = (mats.R(1,1)^2 - r^2*mats.R(2,2)^2)*sin(eps)*cos(eps);
        sigma_d = sqrt(sigma_d);
        sigma_h = sqrt(sigma_h);
        sigma_dh = sqrt(sigma_dh);
        Rk = [sigma_d sigma_dh; sigma_dh sigma_h];
        
        % Covariance update
        P_zz = Rk;
        P_sz = zeros(4,2);
        for j = 1:2*size_x(1) + 1
            P_zz = P_zz + (zeta(:,j) - y_predicted)*(zeta(:,j)-y_predicted)'.*W(j);
            P_sz = P_sz + (sigma_predicted(:,j) - x_hat_predicted)*(zeta(:,j) - y_predicted)'.*W(j);
        end
        K = P_sz/P_zz;
        P_hat_updated = P_hat_updated - K*P_zz*K';
        error_UKF(:,k+1) = y(:,k+1) - y_predicted;
        x_hat_updated(:,k+1) = x_hat_predicted + K*(error_UKF(:,k+1));
        rho_updated(k+1) = max(eig(P_hat_updated));
    end
end