function J = project_f_jacob(x)
    J = [1 0 0 0 0 0 0;...
         0 0 0 -x(7) 0 0 -x(4);...
         0 0 0 1 0 0 0;...
         0 x(7) 0 0 0 0 x(2);...
         0 0 0 0 0 1 0;...
         0 0 0 0 0 0 0;...
         0 0 0 0 0 0 0];
%     J = project_perturb(@project_f, x);
end