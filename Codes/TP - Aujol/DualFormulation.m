function [ub, J, err_u,err_J, niter]=DualFormulation(Image, u0, lambda, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax)
    niter=1;
    
    J=nan(1,itermax);
    err_u=nan(1,itermax);
    err_J=nan(1,itermax);
    
    uold=u0; utilde=u0;
    
    u_x=gradx(utilde); u_y=grady(utilde);
    
    % Initialisation de z
    z1=u_x; z2=u_y;
    normeZ=norm_eps(z1,z2,0);
    nCond=normeZ>1;
    
    z1(nCond)=z1(nCond)./normeZ(nCond);
    z2(nCond)=z2(nCond)./normeZ(nCond);
    
    % Initialiation de I1 et I2
    I1=(Image-c1).^2; I2=(Image-c2).^2;
    
    J(niter)=compute_energy(u0,z1,z2,I1,I2,lambda);  
    % Premiere iteration de l'algorithme
    niter=2; 
    
    z1=z1+tho_z*u_x; z2=z2+tho_z*u_y; 
    normeZ=norm_eps(z1,z2,0);
    nCond=normeZ>1;
    z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
    
    u=uold+tho_u*(div(z1,z2)-lambda*(I1-I2));
    u=min(max(u,0),1);
    utilde=u+theta*(u-uold);
    ub=u>mu;
    
    J(niter)=compute_energy(u,z1,z2,I1,I2,lambda); 
    % Calcul des erreurs
    cond_u=norm(u-uold)/norm(uold);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    err_u(1)=cond_u; err_u(niter)=cond_u;
    err_J(1)=cond_J; err_J(niter)=cond_J;
    %% Boucle while
    while (niter<itermax )%&& cond_u>stop_1 && cond_J>stop_2)
        niter=niter+1;
        uold=u;
        % Dans le cas o� les images sont bruit�es, supposer c1 et c2 connus est faux
        if cinconnu
            c1=sum(sum(Image.*uold))/sum(sum(uold));
            c2=sum(sum(Image.*(1-uold)))/sum(sum(1-uold));

            I1=(Image-c1).^2;
            I2=(Image-c2).^2;
        end
        u_x=gradx(utilde);
        u_y=grady(utilde);
        
        % calcul de z_{k+1}
        z1=z1+tho_z*u_x; z2=z2+tho_z*u_y; 
        normeZ=norm_eps(z1,z2,0);
        nCond=normeZ>1;
        z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
    
        % calcul de u_{k+1}
        u=uold+tho_u*(div(z1,z2)-lambda*(I1-I2));
        u=min(max(u,0),1);
        utilde=u+theta*(u-uold);
        
        ub=u>mu;
        
        J(niter)=compute_energy(u,z1,z2,I1,I2,lambda);
        
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