function main
    
    % Given constant parameters for the model
    params.beta = 40000; % kg/(ms^2) - ballistic coefficient
    params.g = 9.81; % m/s^2 - acceleration due to gravity
    params.Ts = 2; % s - sampling interval
    params.q = 1; % m^2/s^3 - process noise intensity
    params.N = 60; % time steps being considered
    params.dim = 4; % dimension of the state vector
    params.kappa = 1; % UKF Parameter
    
    % Initial values
    init.x0 = 232000; % m 
    init.y0 = 88000; % m
    init.gamma0 = 190 * pi/180; % radians 
    init.v0 = 2290; % m/s
        
    % defining matrices
    mats.Phi = [1 params.Ts 0 0; 0 1 0 0; 0 0 1 params.Ts; 0 0 0 1];
    mats.G = [params.Ts^2/2 0; params.Ts 0; 0 params.Ts^2/2; 0 params.Ts];
    mats.C = [1 0 0 0; 0 0 1 0];
    mats.R = [100^2 0; 0 0.017^2]; % in SI units
    
    % Generating true states
    [params,mats,init,x,y] = true_states(params,mats,init);
%     [mats,y] = measurement(params,mats,x);
%     x_hat_updated = filter_extended_kalman(params,mats,x,y);
    [x_hat_updated_KF,error_KF,rho_predicted_KF,rho_updated_KF] = filter_kalman(init,params,mats,x,y);
    [x_hat_updated_EKF,error_EKF,rho_predicted_EKF,rho_updated_EKF] = filter_extended_kalman(init,params,mats,x,y);
    [x_hat_updated_UKF,error_UKF,rho_predicted_UKF,rho_updated_UKF] = filter_unscented_kalman(init,params,mats,x,y);
    
%     Plots!
    figure
    subplot(2,2,1)
    plot(x(1,:))
    title('True State')
    subplot(2,2,2)
    plot(x(1,:))
    hold on
    plot(x_hat_updated_KF(1,:))
    title('Kalman Filter')
    subplot(2,2,3)
    plot(x(1,:))
    hold on
    plot(x_hat_updated_EKF(1,:))
    title('Extended Kalman Filter')
    legend('Actual Data','EKF')
    subplot(2,2,4)
    plot(x(1,:))
    hold on
    plot(x_hat_updated_UKF(1,:))
    title('Unscented Kalman Filter')
    legend('Actual Data','UKF')

    figure
    subplot(2,2,1)
    plot(x(2,:))
    title('True State')
    subplot(2,2,2)
    plot(x(2,:))
    hold on
    plot(x_hat_updated_KF(2,:))
    title('Kalman Filter')
    subplot(2,2,3)
    plot(x(2,:))
    hold on
    plot(x_hat_updated_EKF(2,:))
    title('Extended Kalman Filter')
    legend('Actual Data','EKF')
    subplot(2,2,4)
    plot(x(2,:))
    hold on
    plot(x_hat_updated_UKF(2,:))
    title('Unscented Kalman Filter')
    legend('Actual Data','UKF')

    figure
    subplot(2,2,1)
    plot(x(3,:))
    title('True State')
    subplot(2,2,2)
    plot(x(3,:))
    hold on
    plot(x_hat_updated_KF(3,:))
    title('Kalman Filter')
    subplot(2,2,3)
    plot(x(3,:))
    hold on
    plot(x_hat_updated_EKF(3,:))
    title('Extended Kalman Filter')
    legend('Actual Data','EKF')
    subplot(2,2,4)
    plot(x(3,:))
    hold on
    plot(x_hat_updated_UKF(3,:))
    title('Unscented Kalman Filter')
    legend('Actual Data','UKF')
    
    figure
    subplot(2,2,1)
    plot(x(4,:))
    title('True State')
    subplot(2,2,2)
    plot(x(4,:))
    hold on
    plot(x_hat_updated_KF(4,:))
    title('Kalman Filter')
    subplot(2,2,3)
    plot(x(4,:))
    hold on
    plot(x_hat_updated_EKF(4,:))
    title('Extended Kalman Filter')
    legend('Actual Data','EKF')
    subplot(2,2,4)
    plot(x(4,:))
    hold on
    plot(x_hat_updated_UKF(4,:))
    title('Unscented Kalman Filter')
    legend('Actual Data','UKF')
    
    figure
    subplot(1,2,1)
    plot(error_KF(1,:))
    hold on
    plot(error_EKF(1,:))
    hold on
    plot(error_UKF(1,:))
    title('Innovation (1)')
    legend('KF','EKF','UKF')
    subplot(1,2,2)
    plot(error_KF(2,:))
    hold on
    plot(error_EKF(2,:))
    hold on
    plot(error_UKF(2,:))
    title('Innovation (2)')
    legend('KF','EKF','UKF')
    
    figure
    plot(sqrt(x(2,:).^2 +x(4,:).^2 ))
    hold on
    plot(sqrt(x_hat_updated_KF(2,:).^2 + x_hat_updated_KF(4,:).^2 ),'r')
    hold on
    plot(sqrt(x_hat_updated_EKF(2,:).^2 + x_hat_updated_EKF(4,:).^2 ),'g')
    hold on
    plot(sqrt(x_hat_updated_UKF(2,:).^2 + x_hat_updated_UKF(4,:).^2 ),'k')
    legend('True','KF','EKF','UKF')
    title('Speed')

    figure
    plot(x(1,:),x(3,:))
    hold on
    plot(x_hat_updated_KF(1,:), x_hat_updated_KF(3,:))
    hold on
    plot(x_hat_updated_EKF(1,:), x_hat_updated_EKF(3,:))
    hold on
    plot(x_hat_updated_UKF(1,:), x_hat_updated_UKF(3,:))
    legend('True','KF','EKF','UKF')
    title('Trajectory')
    
    figure
    subplot(1,2,1)
    plot(rho_predicted_KF)
    hold on
    plot(rho_predicted_EKF)
    hold on
    plot(rho_predicted_UKF)
    title('Spectral radii for the covariances - Predicted')
    legend('KF','EKF','UKF')    
    subplot(1,2,2)
    plot(rho_updated_KF)
    hold on
    plot(rho_updated_EKF)
    hold on
    plot(rho_updated_UKF)
    title('Spectral radii for the covariances - Updated')
    legend('KF','EKF','UKF')

    figure
    subplot(2,2,1)
    plot(x_hat_updated_KF(1,:)-x(1,:))
    hold on
    plot(x_hat_updated_EKF(1,:)-x(1,:))
    hold on
    plot(x_hat_updated_UKF(1,:)-x(1,:))
    title('Error term in KF - State 1')
    legend('KF','EKF','UKF')
    subplot(2,2,2)    
    plot(x_hat_updated_KF(2,:)-x(2,:))
    hold on
    plot(x_hat_updated_EKF(2,:)-x(2,:))
    hold on
    plot(x_hat_updated_UKF(2,:)-x(2,:))
    title('Error term in KF - State 2')
    legend('KF','EKF','UKF')
    subplot(2,2,3)    
    plot(x_hat_updated_KF(3,:)-x(3,:))
    hold on
    plot(x_hat_updated_EKF(3,:)-x(3,:))
    hold on
    plot(x_hat_updated_UKF(3,:)-x(3,:))
    title('Error term in KF - State 3')
    legend('KF','EKF','UKF')
    subplot(2,2,4)    
    plot(x_hat_updated_KF(4,:)-x(4,:))
    hold on
    plot(x_hat_updated_EKF(4,:)-x(4,:))
    hold on
    plot(x_hat_updated_UKF(4,:)-x(4,:))
    title('Error term in KF - State 4')
    legend('KF','EKF','UKF')

% RMSE = sqrt((1/params.N)*sum((x(1,:)-x_hat_updated_KF(1,:)).^2));

end