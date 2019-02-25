function [ub, J, err, niter]=ChanEsedogluNikolova(Image, u0, lambda, mu, tho, sigma, eps, c1, c2, cinconnu, itermax)
    niter=1;
    err=nan(1,itermax);
    J=nan(1,itermax);
    
    uold=u0;
    
    I1=(Image-c1).^2;
    I2=(Image-c2).^2;
    
    ux=gradx(uold);
    uy=grady(uold);
    n_eps=norm_eps(ux,uy,eps);
    
    u=uold+tho*(div(ux./n_eps,uy./n_eps)-lambda*(I1-I2));
    u=min(max(u,0),1);
    ub=u>mu;
    
    err(niter)=norm(u-uold);
    J(niter)=compute_energy_smooth(ub,I1,I2,lambda,eps); 
    
    while (niter<itermax && norm(u-uold)>sigma)
        niter=niter+1;
        uold=u;

        if cinconnu
            c1=sum(sum(Image.*uold))/sum(sum(uold));
            c2=sum(sum(Image.*(1-uold)))/sum(sum(1-uold));

            I1=(Image-c1).^2;
            I2=(Image-c2).^2;
        end
    
        ux=gradx(u);
        uy=grady(u);
        n_eps=norm_eps(ux,uy,eps);
    
        u=uold+tho*(div(ux./n_eps,uy./n_eps)-lambda*(I1-I2));
        u=min(max(u,0),1);
        ub=u>mu;
        
        err(niter)=norm(u-uold);
        J(niter)=compute_energy_smooth(ub,I1,I2,lambda,eps); 
    end
    err(isnan(err))=[];
    J(isnan(J))=[];
end