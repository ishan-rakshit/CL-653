function main
    % Generating true states
    [params,mats,x] = true_states;
    [mats,y] = measurement(params,mats,x);
    x_hat_updated = filter_extended_kalman(params,mats,x,y);
    
end