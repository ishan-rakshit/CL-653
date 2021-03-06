function [c1,c2,rho] = return_rho(x_3)
    if x_3 < 9144
        c1 = 1.227;
        c2 = 1.093e-04;
    else
        c1 = 1.754;
        c2 = 1.49e-04;
    end
    rho = c1*exp(-1*c2*x_3);
end
