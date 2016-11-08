function unscented_transform(x,n_x,x_bar,P_x,y)
% calculating the first two moments of y using the unscented transform (UT)
    chi = zeros(length(x(:,1)),2*n_x+1); % sigma points
    W = zeros(2*n_x+1,1); % weights of the sigma points
    chi(:,1) = x_bar;
    W(1,1) = kappa/(n_x + kappa);
    y(1,1) = return_g(chi(:,1));
    for i = 2:n_x
        chi(:,i) = x_bar + nkp_mat(:,i);
        W(i,1) = 1/(2*(n_x + 1));
        y(i,1) = return_g(chi(:,i)); % propagation of points through non-linear function
    end
    for i = n_x + 1:2*n_x
        chi(:,i) = x_bar - nkp_mat(:,i);
        W(i,1) = 1/(2*(n_x + 1));
        y(i,1) = return_g(chi(:,i)); % propagation of points through non-linear function
    end
    y_bar = sum(W.*y);
    P_bar = sum(W.*(y-y_bar).*transpose(y-y_bar));
    
    % Implementation of UKF
end