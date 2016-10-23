function project_main
    randn('state', 0);
    sig_1 = 0.2;
    sig_2 = 7E-3;
    sig_r = 50;
    sig_theta = 0.1;
    sig_phi = 0.1;
    T_s = 0.1;
    N = 1000;
    
    Q_d = diag([0 sig_1^2 0 sig_1^2 0 sig_1^2 sig_2^2]);
    w_d = sqrt(Q_d)*randn(7, N);
    disp('Q_d (sample) is:');
    disp(cov(w_d'));
    R = diag([sig_r^2 sig_theta^2 sig_phi^2]);
    v = sqrt(R)*randn(3, N);
    disp('R (sample) is:');
    disp(cov(v'));
    
    x_0 = [1000; 0; 2650; 150; 200; 0; 3];
    y = [];
    for k = 1:N
       [t, x_temp] = ode45(@project_f, [0, T_s], x_0);
       x(:,k) =  x_temp(end,:)' + w_d(:,k);
       t = 0; % redundant
       y(:,k) = project_g(t, x(:,k)) + v(:,k);
       x_0 = x(:,k);
    end
    
    % KALMAN FILTER
    x_hat_0_KF = zeros(7,1);
    P_KF = eye(7);
    B_d = eye(7);
    x_hat_KF = [];
    e_KF = [];
    tr_KF = [];
    for k = 1:N
        % Local linearisation
        A = project_perturb(@project_f, x_hat_0_KF); % x(k-1|k-1)
        Phi = expm(A*T_s);
        % A is not coming out to be invertible; can't apply this
        % Gam_d = (A/(Phi-eye(7)))*B_d;
        Gam_d = integral(@(t) expm(A*(T_s-t))*B_d,0,T_s,'ArrayValued',true);
        C = project_perturb(@project_g, x_hat_0_KF);
        % Predict
        x_hat_0_KF = Phi*x_hat_0_KF;
        P_KF = Phi*P_KF*Phi'+Gam_d*Q_d*Gam_d';
        L = (P_KF*C')*inv(C*P_KF*C'+R);
        % Estimate
        e_KF(:,k) = (y(:,k)-C*x_hat_0_KF);
        x_hat_0_KF = x_hat_0_KF + L*e_KF(:,k); % x(k|k)
        x_hat_KF(:,k) = x_hat_0_KF;
        P_KF = (eye(7)-L*C)*P_KF; % P(k|k)
        tr_KF(k) = trace(P_KF);
    end
    
    % EXTENDED KALMAN FILTER
    x_hat_0_EKF = zeros(7,1);
    P_EKF = eye(7);
    B_d = eye(7);
    x_hat_EKF = [];
    e_EKF = [];
    tr_EKF = [];
    for k = 1:N
        % Local linearisation
        A = project_perturb(@project_f, x_hat_0_EKF); % x(k-1|k-1)
        Phi = expm(A*T_s);
        Gam_d = integral(@(t) expm(A*(T_s-t))*B_d,0,T_s,'ArrayValued',true);
        C = project_perturb(@project_g, x_hat_0_EKF);
        % Predict
        [t, x_temp] = ode45(@project_f, [0, T_s], x_hat_0_EKF);
        x_hat_EKF(:,k) =  x_temp(end,:)';
        P_EKF = Phi*P_EKF*Phi'+Gam_d*Q_d*Gam_d';
        L = (P_EKF*C')*inv(C*P_EKF*C'+R);
        % Estimate
        e_EKF(:,k) = (y(:,k)-C*x_hat_0_EKF);
        x_hat_0_EKF = x_hat_0_EKF + L*e_EKF(:,k); % x(k|k)
        x_hat_EKF(:,k) = x_hat_0_EKF;
        P_EKF = (eye(7)-L*C)*P_EKF; % P(k|k)
        tr_EKF(k) = trace(P_EKF);
    end
    
    % PLOT
    x_0 = [1000; 0; 2650; 150; 200; 0; 3];
    x_hat_0_KF = zeros(7,1);
    figure;
    % state 1
    subplot(4,4,1)
    plot(0:N, [x_0(1), x(1,:)], 'r', 0:N, [x_hat_0_KF(1), x_hat_KF(1,:)], 'g', 0:N, [x_hat_0_EKF(1), x_hat_EKF(1,:)], 'b')
    title('x_1')
    legend('True', 'KF', 'EKF')
    % state 2
    subplot(4,4,2)
    plot(0:N, [x_0(2), x(2,:)], 'r', 0:N, [x_hat_0_KF(2), x_hat_KF(2,:)], 'g', 0:N, [x_hat_0_EKF(2), x_hat_EKF(2,:)], 'b')
    title('x_2')
    % state 3
    subplot(4,4,3)
    plot(0:N, [x_0(3), x(3,:)], 'r', 0:N, [x_hat_0_KF(3), x_hat_KF(3,:)], 'g', 0:N, [x_hat_0_EKF(3), x_hat_EKF(3,:)], 'b')
    title('x_3')
    % state 4
    subplot(4,4,4)
    plot(0:N, [x_0(4), x(4,:)], 'r', 0:N, [x_hat_0_KF(4), x_hat_KF(4,:)], 'g', 0:N, [x_hat_0_EKF(4), x_hat_EKF(4,:)], 'b')
    title('x_4')
    % state 5
    subplot(4,4,5)
    plot(0:N, [x_0(5), x(5,:)], 'r', 0:N, [x_hat_0_KF(5), x_hat_KF(5,:)], 'g', 0:N, [x_hat_0_EKF(5), x_hat_EKF(5,:)], 'b')
    title('x_5')
    % state 6
    subplot(4,4,6)
    plot(0:N, [x_0(6), x(6,:)], 'r', 0:N, [x_hat_0_KF(6), x_hat_KF(6,:)], 'g', 0:N, [x_hat_0_EKF(6), x_hat_EKF(6,:)], 'b')
    title('x_6')
    % state 7
    subplot(4,4,7)
    plot(0:N, [x_0(7), x(7,:)], 'r', 0:N, [x_hat_0_KF(7), x_hat_KF(7,:)], 'g', 0:N, [x_hat_0_EKF(7), x_hat_EKF(7,:)], 'b')
    title('x_7')
    % output 1
    subplot(4,4,8)
    plot(1:N, y(1,:), 'r')
    title('y_1')
    % output 2
    subplot(4,4,9)
    plot(1:N, y(2,:), 'r')
    title('y_2')
    % output 3
    subplot(4,4,10)
    plot(1:N, y(3,:), 'r')
    title('y_3')
    % innovation 1
    subplot(4,4,11)
    plot(1:N, e_KF(1,:), 'g', 1:N, e_EKF(1,:), 'b')
    title('e_1')
    % innovation 2
    subplot(4,4,12)
    plot(1:N, e_KF(2,:), 'g', 1:N, e_EKF(2,:), 'b')
    title('e_2')
    % innovation 3
    subplot(4,4,13)
    plot(1:N, e_KF(3,:), 'g', 1:N, e_EKF(3,:), 'b')
    title('e_3')
    % Trace
    subplot(4,4,14)
    plot(1:N, tr_KF, 'g', 1:N, tr_EKF, 'b')
end