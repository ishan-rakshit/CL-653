function f = return_f(xk,params)
    temp_vec = [xk(2) xk(4)]';
    f = -0.5*(params.g*return_rho(xk(3))/params.beta)*norm(temp_vec)*temp_vec;
end