function J = project_g_jacob(x)
    J = project_perturb(@project_g, x);
end