function out = return_F_XUD(x, k, mats, params, w)
    out = return_psi(x(:,k),mats,params) + mats.G*[0 -params.g]' + w(:,k);
end