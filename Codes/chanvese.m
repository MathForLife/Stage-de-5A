function [phib, J, step, err, c1,c2, niter]=chanvese(Image, mask, lambda, sigma, eps, eta, n, order, itermax)
    niter=1;
    J=nan(1,itermax);
    step=nan(1,itermax);
    err=nan(1,itermax);
    
    phi=signed_distance_from_mask(mask);
    
    Hold=Heavyside_eta(phi, eta, order);
    dold=delta_eta(phi,eta, order);
    
    c1=sum(sum(Image.*Hold))/sum(sum(Hold));
    c2=sum(sum(Image.*(1-Hold)))/sum(sum(1-Hold));
    
    I1=(Image-c1).^2;
    I2=(Image-c2).^2;
    
    phix=gradx(phi);
    phiy=grady(phi);
    n_eps=norm_eps(phix,phiy,eps);
    
    GradJ=dold.*(div(phix./n_eps,phiy./n_eps)-lambda*(I1-I2));
    tho=0.5/max(max(GradJ));
    
    phi=phi+tho*GradJ;
    phib=phi>0;
    
    Hn=Heavyside_eta(phi, eta, order);
    dn=delta_eta(phi,eta, order);

    step(niter)=tho;
    err(niter)=norm(Hn-Hold);
    J(niter)=compute_energy_smooth(phib,I1,I2,lambda,eps); 
    
    while (niter<itermax && norm(Hn-Hold)>sigma)
        
        niter=niter+1;
        
        Hold=Hn;
        dold=dn;
        
        c1=sum(sum(Image.*Hold))/sum(sum(Hold));
        c2=sum(sum(Image.*(1-Hold)))/sum(sum(1-Hold));
        
        I1=(Image-c1).^2;
        I2=(Image-c2).^2;
        
        phix=gradx(phi);
        phiy=grady(phi);
        n_eps=norm_eps(phix,phiy,eps);
    
        GradJ=dold.*(div(phix./n_eps,phiy./n_eps)-lambda*(I1-I2));
        
        tho=0.5/max(max(GradJ));
        phi=phi+tho*GradJ;
        if mod(niter,n)==0
            phi=signed_distance_from_mask(phi>0);
        end
        phib=phi>0;
        
        Hn=Heavyside_eta(phi, eta, order);
        dn=delta_eta(phi,eta, order);
   
        step(niter)=tho;
        err(niter)=norm(Hn-Hold);
        J(niter)=compute_energy_smooth(phib,I1,I2,lambda,eps);
    end
    step(isnan(step))=[];
    err(isnan(err))=[];
    J(isnan(J))=[];
end