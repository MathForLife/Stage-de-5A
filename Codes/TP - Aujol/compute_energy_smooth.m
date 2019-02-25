function Energy=compute_energy_smooth(u,I1,I2,lambda,eps)
    ux=gradx(u);
    uy=grady(u);
    n_eps=norm_eps(ux,uy,eps);
    
    Energy=sum(sum(n_eps+lambda*I1.*u+lambda*I2.*(1-u)));
end