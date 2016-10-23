function rho = return_rho(x_3)
    if x_3 < 9144
        rho = 1.227*exp(-1.093e-04 * x_3);
    else
        rho = 1.754*exp(-1.49e-04 * x_3);
    end
end
