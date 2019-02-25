function n_eps=norm_eps(Ix,Iy,eps)
    n_eps=sqrt(Ix.^2+Iy.^2+eps.^2);
end