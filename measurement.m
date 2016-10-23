function y = measurement(params,x)
    [~,N] = size(x);
    
    % measurement noise
    R = [100^2 0; 0 0.017^2]; % in SI units
    v = mvnrnd(zeros(1,2),R,params.N)';
    y = zeros(2,params.N);
   
    for k = 1:N
        temp_c = x(1,k) + 1i * x(3,k);
%         x3 = x(3,k);
        y(:,k) = [abs(temp_c), angle(temp_c)]' + v(:,k);
    end
    
    
    
%     % checking v
%     temp_v = 100*randn(2,params.N);
%     temp_y = zeros(2,params.N);
%     for k = 1:N
%         temp_c = x(1,k) + 1i * x(3,k);
% %         x3 = x(3,k);
%         temp_y(:,k) = [abs(temp_c), angle(temp_c)]' + temp_v(:,k);
%     end

end