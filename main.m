function main
    % Generating true states
    [params,mats,x] = true_states;
    [mats,y] = measurement(params,mats,x);
%     x_hat_updated = filter_extended_kalman(params,mats,x,y);
    x_hat_updated = filter_kalman(params,mats,x,y);
    
    
    % Let the plotting begin!
    figure
    subplot(4,1,1)
    plot(x(1,:))
    hold on
    plot(x_hat_updated(1,:))
%     hold on
    title('x1 comparison')
    subplot(4,1,2)
    plot(x(2,:),'b')
    hold on
    plot(x_hat_updated(2,:),'r--')
%     hold on
    title('x2 comparison')
    subplot(4,1,3)
    plot(x(3,:),'b')
    hold on
    plot(x_hat_updated(3,:),'r--')
%     hold on
    title('x3 comparison')
    subplot(4,1,4)
    plot(x(4,:),'b')
    hold on
    plot(x_hat_updated(4,:),'r--')
%     hold on
    title('x4 comparison')
    
end