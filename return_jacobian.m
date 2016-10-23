function Psi_k = return_jacobian(params, x_vector)
%     x1 = x_vector(1);
    x2 = x_vector(2);
    x3 = x_vector(3);
    x4 = x_vector(4);
    [~,c2,rho] = return_rho(x3);
    F_k(1,1) = 0; 
    F_k(1,2) = -0.5*(params.g/params.beta)*rho*((2*(x2)^2 + (x4)^2)/(sqrt((x2)^2 + (x4)^2))); 
    F_k(1,3) = 0.5*(params.g/params.beta)*c2*rho*((x2)*(sqrt((x2)^2 + (x4)^2)));  
    F_k(1,4) = -0.5*(params.g/params.beta)*rho*((x2*x4)/(sqrt((x2)^2 + (x4)^2))); 
    F_k(2,1) = 0; 
    F_k(2,2) = F_k(1,4); 
    F_k(2,3) = 0.5*(params.g/params.beta)*c2*rho*((x4)*(sqrt((x2)^2 + (x4)^2))); 
    F_k(2,4) = -0.5*(params.g/params.beta)*rho*(((x2)^2 + 2*(x4)^2)/(sqrt((x2)^2 + (x4)^2))); 
    
    % computing the jacobian Psi_k 
    Psi_k = mats.Phi + mats.G*F_k;
end