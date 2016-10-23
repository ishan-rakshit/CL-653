function Psi = return_jacobian(X_star, U_star, D_star, params, means)
% Input: name of function created above, X?,U ?,D? (Here X?; U?; D? are the values of X, U and D at which the Jacobian matrices have to be evaluated.) 
% Output: three Jacobian matrices, namely A; B; H. 
    % initializing Jacobian matrices
    size_of_vector = 2;
    A = zeros(size_of_vector,size_of_vector);
%     B = zeros(size_of_vector,1);
%     H = zeros(size_of_vector,1);
    
    % perturbation values
    delta_X = 0.001;
    delta_X1 = [delta_X; 0];
    delta_X2 = [0; delta_X];
    delta_U = 0.001;
    delta_D = 0.001;
    
    % creating matrices
    A(:,1) = (0.5/delta_X)*(ode_rhs(X_star + delta_X1, U_star, D_star, params, means) - ode_rhs(X_star - delta_X1, U_star, D_star, params, means));
    A(:,2) = (0.5/delta_X)*(ode_rhs(X_star + delta_X2, U_star, D_star, params, means) - ode_rhs(X_star - delta_X2, U_star, D_star, params, means));
    B = (0.5/delta_U)*(ode_rhs(X_star, U_star + delta_U, D_star, params, means) - ode_rhs(X_star, U_star - delta_U, D_star, params, means));
    H = (0.5/delta_D)*(ode_rhs(X_star, U_star, D_star + delta_D, params, means) - ode_rhs(X_star, U_star, D_star - delta_D, params, means));
end