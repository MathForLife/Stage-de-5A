function [ub, J, err_u,err_J, niter]=ChanEsedogluNikolova(Image, u0, lambda, mu, tho, eps, stop_1, stop_2, c1, c2, cinconnu, itermax)
    niter=1;
    
    err_u=nan(1,itermax);
    err_J=nan(1,itermax);
    J=nan(1,itermax);
    
    uold=u0;
    
    I1=(Image-c1).^2;
    I2=(Image-c2).^2;
    
    J(niter)=compute_energy_smooth(u0,I1,I2,lambda,eps);
    %% Première itération
    niter=2;
    
    ux=gradx(uold);
    uy=grady(uold);
    n_eps=norm_eps(ux,uy,eps);
    
    u=uold+tho*(div(ux./n_eps,uy./n_eps)-lambda*(I1-I2));
    u=min(max(u,0),1);
    ub=u>mu;
    
    J(niter)=compute_energy_smooth(ub,I1,I2,lambda,eps); 
    
    % Calcul des erreurs
    cond_u=norm(u-uold)/norm(uold);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    err_u(1)=10; err_u(niter)=cond_u;
    err_J(1)=10; err_J(niter)=cond_J;
    
    while (niter<itermax && cond_u>stop_1 && cond_J>stop_2)
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
        
        J(niter)=compute_energy_smooth(ub,I1,I2,lambda,eps);
        
        % Calcul des erreurs
        cond_u=norm(u-uold)/norm(uold);
        cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
        if J(niter)>J(niter-1)
            niter=niter-1;
            
            u=uold;
            tho=tho/2;
            fprintf('\tho= %5.3f, niter=%d\n',tho,niter)
        else
            err_u(niter)=cond_u;
            err_J(niter)=cond_J;
        end
        
    end
    err_u(isnan(err_u))=[];
    err_J(isnan(err_J))=[];
    J(isnan(J))=[];
end