% Submitted by Ankita Humne (120010012)
% Project CL653

clc 
close all
clear all

T = 2;%sec
g = 9.81; %m/s2
beta = 40000; %kg/ms2
q = 1; %m2/s3
N = 60;
randn('state',0)

phi = [1 T 0 0 ; 0 1 0 0; 0 0 1 T; 0 0 0 1];
G = [T^2/2 0; T 0; 0 T^2/2; 0 T];
Q = q*[T^3/3 T^2/2 0 0; T^2/2 T 0 0;0 0 T^3/3 T^2/2; 0 0 T^2/2 T];
R = [100^2 0; 0 0.017^2];
C = [1 0 0 0; 0 0 1 0];

w = mvnrnd([0 0 0 0],Q,N); %process noise
% v = mvnrnd([0 0],R,N); %measurement noise
x = zeros(4,N);
y = zeros(2,N);
x(:,1) = [232000; -2.2552e3; 88000; -397.6543] ;
rho = 1.754*exp(-1.49e-4*x(3,1));
v = mvnrnd([0 0],R);
r = sqrt(x(3,1)^2 + x(1,1)^2) + v(1); 

e_KF = 190*pi/180;%atan2(x(3,1),x(1,1)) + v(2);
R_c(1,1) = R(1,1)^2*(cos(e_KF))^2 + r^2*R(2,2)*(sin(e_KF))^2;
R_c(2,2) = R(1,1)^2*(sin(e_KF))^2 + r^2*R(2,2)*(cos(e_KF))^2;
R_c(1,2) = (R(1,1)^2 - r^2*R(2,2)^2)*sin(e_KF)*cos(e_KF);
R_c(2,1) = R_c(1,2);
v = mvnrnd([0 0],R_c);
    
%% True State Generation    
for i = 1:N-1
    if x(3,i) < 9144
        rho = 1.227*exp(-1.093e-4*x(3,i));
    elseif x(3,i) >= 9144
        rho = 1.754*exp(-1.49e-4*x(3,i));
    end
    f = -[x(2,i);x(4,i)].*0.5*(g*rho/beta)*sqrt((x(2,i))^2 + (x(4,i))^2);
    x(:,i+1) = phi*x(:,i) + G*f + G*[0; -g] + w(i,:)';
    
    y(:,i) = C*x(:,i) + v'; 
    
    r = sqrt(y(1,i)^2 + y(2,i)^2); 
    e_KF = atan2(y(2,i),y(1,i));
    R_c(1,1) = R(1,1)^2*(cos(e_KF))^2 + r^2*R(2,2)*(sin(e_KF))^2;
    R_c(2,2) = R(1,1)^2*(sin(e_KF))^2 + r^2*R(2,2)*(cos(e_KF))^2;
    R_c(1,2) = (R(1,1)^2 - r^2*R(2,2)^2)*sin(e_KF)*cos(e_KF);
    R_c(2,1) = R_c(1,2);
    v = mvnrnd([0 0],R_c);
end

    v = mvnrnd([0 0],R_c);

y(:,N) = C*x(:,N)+v';


r = sqrt(y(1,1)^2 + y(2,1)^2); 
el = 190*pi/180;%atan2(y(2,1),y(1,1));
sd = R(1,1)^2*(cos(el))^2 + r^2*R(2,2)*(sin(el))^2;
sh = R(1,1)^2*(sin(el))^2 + r^2*R(2,2)*(cos(el))^2;
sdh = (R(1,1)^2 - r^2*R(2,2)^2)*sin(el)*cos(el)  ; 
P_init = [sd -sd/T sdh -sdh/T; -sd/T 2*sd/T^2 -sdh/T 2*sdh/T^2; sdh -sdh/T sh -sh/T; -sdh/T 2*sdh/T^2 -sh/T 2*sh/T^2]; %covariance matrix

%% Kalman filter implementation
x_eq = [-295000; 0; -5000; 0 ];%x(:,1);
dx_KF = zeros(4,N);
y_eq = C*x_eq; 
F = f_jacobian(x(:,1));
x_KF = zeros(4,N);
x_KF(:,1) = x(:,1);
dx_KF(:,1) = x_KF(:,1) - x_eq;
e_KF = zeros(2,N-1);
rho_pred = zeros(1,N);
rho_upd = zeros(1,N);
P = P_init;
rho_pred(1) = max(eig(P));
rho_upd(1) = max(eig(P));
f_eq = -[x_eq(2);x_eq(4)].*0.5*(g*rho/beta)*sqrt((x_eq(2))^2 + (x_eq(4))^2);

for i = 1:N-1
    
    
    r = sqrt(y(2,i)^2 + y(1,i)^2); 
    eta_KF = atan2(y(2,i),y(1,i));
    R_c(1,1) = R(1,1)^2*(cos(eta_KF))^2 + r^2*R(2,2)*(sin(eta_KF))^2;
    R_c(2,2) = R(1,1)^2*(sin(eta_KF))^2 + r^2*R(2,2)*(cos(eta_KF))^2;
    R_c(1,2) = (R(1,1)^2 - r^2*R(2,2)^2)*sin(eta_KF)*cos(eta_KF);
    R_c(2,1) = R_c(1,2);
    
%     f = -[x_KF(2,i);x_KF(4,i)].*0.5*(g*rho/beta)*sqrt((x_KF(2,i))^2 + (x_KF(4,i))^2);
    
%     x_pred = phi*(x_KF(:,i)-x_eq) + G*(f-f_eq) + G*[0; -g] ;
    x_pred = (phi + G*F)*dx_KF(:,i)+ G*[0;-g];
    P = (phi + G*F)*P*(phi + G*F)' + Q;
    rho_pred(i+1) = max(eig(P));
    L = P*C'*inv(C*P*C'+R);
    e_KF(:,i) = y(:,i)- C*x_pred - y_eq;
    dx_KF(:,i+1) = x_pred + L*e_KF(:,i);
    x_KF(:,i+1) = x_eq + dx_KF(:,i+1);
    P = (eye(4)-L*C)*P;
    rho_upd(i+1) = max(eig(P));
    
end



%% Extended kalman filter implementation
x_EKF = zeros(4,N);
x_EKF(:,1) = x(:,1);
e_EKF = zeros(2,N-1);
P = P_init;

for i = 1:N-1
    if x_EKF(3,i) < 9144
        rho = 1.227*exp(-1.093e-4*x_EKF(3,i));
    elseif x_EKF(3,i) >= 9144
        rho = 1.754*exp(-1.49e-4*x_EKF(3,i));
    end
    F = f_jacobian(x_EKF(:,i));
    f = -[x_EKF(2,i);x_EKF(4,i)].*0.5*(g*rho/beta)*sqrt((x_EKF(2,i))^2 + (x_EKF(4,i))^2);
    x_pred = phi*x_EKF(:,i) + G*f + G*[0; -g] ;
%     x_pred = (phi + G*F)*x_EKF(:,i) + G*[0;-g];
    P = (phi + G*F)*P*(phi + G*F)' + Q;
    
    r = sqrt(y(2,i)^2 + y(1,i)^2); 
    eta_EKF = atan2(y(2,i),y(1,i));
    R_c(1,1) = R(1,1)^2*(cos(eta_EKF))^2 + r^2*R(2,2)*(sin(eta_EKF))^2;
    R_c(2,2) = R(1,1)^2*(sin(eta_EKF))^2 + r^2*R(2,2)*(cos(eta_EKF))^2;
    R_c(1,2) = (R(1,1)^2 - r^2*R(2,2)^2)*sin(eta_EKF)*cos(eta_EKF);
    R_c(2,1) = R_c(1,2);
    
    L = P*C'/(C*P*C' + R_c);
    e_EKF(:,i) = y(:,i)-C*x_pred; 
    x_EKF(:,i+1) = x_pred + L*e_EKF(:,i);
    P = (eye(4) - L*C)*P;

end


%% Unscented Kalman Filter
s = size(x);
nx = s(1);
% x_bar = [mean(x(1,:));mean(x(2,:));mean(x(3,:));mean(x(4,:))];

kappa = 1;
W = ones(1,2*nx+1).*1/(2*(nx + 1));
W(1) = kappa/(nx+kappa);
x_UKF = zeros(4,N);

 

% sigma(:,1)= x_bar;
P = P_init;
Ps = P.*(nx + kappa);
Psq = sqrtm(Ps);
% for i=2:nx+1
%     sigma(:,i) = x_bar + Psq(:,(i-1));
% end
% 
% for i = nx+2 : 2*nx+1
%     sigma(:,i) = x_bar - Psq(:,(i-nx-1));
% end
% x1 = zeros(4,1);
% for i = 1:2*nx+1
%     x1 = x1 + sigma(:,i).*W(i); 
% end
x_UKF(:,1) = x(:,1);
zeta = zeros(2,2*nx+1);
for i = 1: N-1
%     F = f_jacobian(sigma(:,i));
    sigma = zeros(4,2*nx+1);
    sigma(:,1) = x_UKF(:,i);
    for j=2:nx+1
        sigma(:,j) = x_UKF(:,i) + Psq(:,(j-1));
    end
    for j = nx+2 : 2*nx+1
        sigma(:,j) = x_UKF(:,i) - Psq(:,(j-nx-1));
    end
    sigma_pred = zeros(4,2*nx+1);
    x_pred = zeros(4,1);
    for j = 1: 2*nx +1
        f = -[sigma(2,j);sigma(4,j)].*0.5*(g*rho/beta)*sqrt((sigma(2,j))^2 + (sigma(2,j))^2);
        sigma_pred(:,j) = phi*sigma(:,j) + G*f + G*[0; -g] ;
%         sigma_pred(:,j) = (phi + G*F)*sigma(:,j) + G*[0;-g];
        x_pred = x_pred + sigma_pred(:,j).*W(j);
    end
    P = Q;
    for j = 1:2*nx+1
        P = P+(sigma_pred(:,j)-x_pred)*(sigma_pred(:,j)-x_pred)'.*W(j);
    end
    for j = 1:2*nx +1
        zeta(:,j) = C*sigma_pred(:,j);
    end
    y_pred = zeros(2,1);
    for j = 1:2*nx+1
        y_pred = y_pred + zeta(:,j).*W(j);
    end
    
    
    r = sqrt(y(2,i)^2 + y(1,i)^2); 
    eta_EKF = atan2(y(2,i),y(1,i));
    R_c(1,1) = R(1,1)^2*(cos(eta_EKF))^2 + r^2*R(2,2)*(sin(eta_EKF))^2;
    R_c(2,2) = R(1,1)^2*(sin(eta_EKF))^2 + r^2*R(2,2)*(cos(eta_EKF))^2;
    R_c(1,2) = (R(1,1)^2 - r^2*R(2,2)^2)*sin(eta_EKF)*cos(eta_EKF);
    R_c(2,1) = R_c(1,2);
    
    Pzz = R_c;
    Psz = zeros(4,2);
    for j = 1: 2*nx+1
        Pzz  = Pzz + (zeta(:,j)-y_pred)*(zeta(:,j)-y_pred)'.*W(j);
        Psz = Psz + (sigma_pred(:,j)-x_pred)*(zeta(:,j)-y_pred)'.*W(j);
    end
    
    K = Psz/(Pzz);
    x_UKF(:,i+1) = x_pred + K*(y(:,i+1) - y_pred);
    P = P - K*Pzz*K';
    
     
end
    
    
    
    

figure
plot(sqrt(x(2,:).^2 +x(4,:).^2 ))
hold on
plot(sqrt(x_KF(2,:).^2 +x_KF(4,:).^2 ),'r')
hold on
plot(sqrt(x_EKF(2,:).^2 +x_EKF(4,:).^2 ),'g')
hold on
plot(sqrt(x_UKF(2,:).^2 +x_UKF(4,:).^2 ),'k')
legend('True','KF','EKF','UKF')
title('Speed')


figure
plot(x(1,:)./1000,x(3,:)./1000)
hold on
plot(x_KF(1,:)./1000,x_KF(3,:)./1000,'r')
hold on
plot(x_EKF(1,:)./1000,x_EKF(3,:)./1000,'g')
hold on
plot(x_UKF(1,:)./1000,x_UKF(3,:)./1000,'k')
legend('True','KF','EKF','UKF')
title('Trajectory')
% 
% figure 
% plot(x(2,:))
% hold on
% plot(x_KF(2,:),'r')
% hold on
% plot(x_EKF(2,:),'g')
% legend('True','KF','EKF')
% title('v_x')
% 
% figure
% plot(x(4,:))
% hold on
% plot(x_KF(4,:),'r')
% hold on
% plot(x_EKF(4,:),'g')
% legend('True','KF','EKF')
% title('v_y')

figure
plot(e_KF(1,:))
hold on
plot(e_EKF(1,:),'r')
legend('KF','EKF')
title('Innovation_y')


figure
plot(e_KF(2,:))
hold on
plot(e_EKF(2,:),'r')
legend('KF','EKF')
title('Innovation_y')
% 



% figure
% plot(rho_pred)
% hold on
% plot(rho_upd,'r')
% legend('Predicted','Updated')
% title('Spectral radius of Covariance matrix')