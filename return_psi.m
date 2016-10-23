function psi = return_psi(xk,mats,params)
    psi = mats.Phi*xk + mats.G*return_f(xk,params);
end