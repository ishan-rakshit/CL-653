function [ F ] = f_jacobian( x )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
T = 2;%sec
g = 9.81; %m/s2
beta = 40000; %kg/ms2

if x(3) < 9144
    rho = 1.227*exp(-1.093e-4*x(3));
    c_2 = 1.093e-4;
elseif x(3) >= 9144
    rho = 1.754*exp(-1.49e-4*x(3));
    c_2 = 1.49e-4;
end

F = zeros(2,4);
F(1,2) = -0.5*(g*rho/beta)*(2*x(2)^2+x(4)^2)/sqrt(x(2)^2+x(4)^2);
F(1,3) = 0.5*(c_2*g*rho/beta)*x(2)*sqrt(x(2)^2+x(4)^2);
F(1,4) = -0.5*(g*rho/beta)*(x(2)*x(4))/sqrt(x(2)^2+x(4)^2);
F(2,2) = F(1,4);
F(2,3) = 0.5*(c_2*g*rho/beta)*x(4)*sqrt(x(2)^2+x(4)^2);
F(2,4) = -0.5*(g*rho/beta)*(2*x(4)^2+x(2)^2)/sqrt(x(2)^2+x(4)^2);

end

