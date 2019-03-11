function [phib, J, err_u,err_J, niter]=ChanVese(Image, mask, lambda, eps, eta,stop_1,stop_2,c1,c2, n, order, itermax)
    niter=1;
    J=nan(1,itermax);
    err_u=nan(1,itermax);
    err_J=nan(1,itermax);
    
    phi=signed_distance_from_mask(mask);
    
    Hold=Heavyside_eta(phi, eta, order);
    dold=delta_eta(phi,eta, order);
    
    I1=(Image-c1).^2;
    I2=(Image-c2).^2;
    
    J(niter)=compute_energy_smooth(mask,I1,I2,lambda,eps); 
    
    %% Première itération
    niter=2;
    phix=gradx(phi);
    phiy=grady(phi);
    n_eps=norm_eps(phix,phiy,eps);
    
    GradJ=dold.*(div(phix./n_eps,phiy./n_eps)-lambda*(I1-I2));
    tho=0.5/max(max(GradJ));
    
    phi=phi+tho*GradJ;
    phib=phi>0;
    
    Hn=Heavyside_eta(phi, eta, order);
    dn=delta_eta(phi,eta, order);

    J(niter)=compute_energy_smooth(phib,I1,I2,lambda,eps); 
    % Calcul des erreurs
    cond_u=norm(Hn-Hold)/norm(Hold);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    err_u(1)=10; err_u(niter)=cond_u;
    err_J(1)=10; err_J(niter)=cond_J;
    while (niter<itermax && cond_u>stop_1 && cond_J>stop_2)
        
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
   
        J(niter)=compute_energy_smooth(phib,I1,I2,lambda,eps);
        
        % Calcul des erreurs
        cond_u=norm(u-uold)/norm(uold);
        cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));

        err_u(niter)=cond_u;
        err_J(niter)=cond_J;
    end
    err_u(isnan(err_u))=[];
    err_J(isnan(err_J))=[];
    J(isnan(J))=[];
end