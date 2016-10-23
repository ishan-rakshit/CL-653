function J = project_perturb(project_dot, x)
    epsi = 0.01*eye(size(x,1)); % each column in your perturbation along one direction
    J = [];
    t = 0; % redundant
    for k = 1:size(x,1)
        J(:,k) = (project_dot(t, x+epsi(:,k))-project_dot(t, x-epsi(:,k)))/(2*norm(epsi(:,k)));
    end
end