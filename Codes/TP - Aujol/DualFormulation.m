function [ub, J1, J2, err, niter]=DualFormulation(Image, u0, lambda, mu, tho, theta, sigma,stopCond, c1, c2, cinconnu, itermax)
    niter=1;
    J1=nan(1,itermax);
    J2=nan(1,itermax);
    err=nan(1,itermax);
    
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
    
    % Premiere iteration de l'algorithme
    z1=z1+sigma*u_x; z2=z2+sigma*u_y; 
    normeZ=norm_eps(z1,z2,0);
    nCond=normeZ>1;
    z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
    
    u=uold+tho*(div(z1,z2)-lambda*(I1-I2));
    u=min(max(u,0),1);
    utilde=u+theta*(u-uold);
    ub=u>mu;
    
    err(niter)=norm(u-uold);
    J1(niter)=compute_energy_smooth(ub,I1,I2,lambda,0); 
    J2(niter)=compute_energy(u,z1,z2,I1,I2,lambda); 
    
    %% Boucle while
    while (niter<itermax )
        niter=niter+1;
        uold=u;
        % Dans le cas où les images sont bruitées, supposer c1 et c2 connus est faux
        if cinconnu
            c1=sum(sum(Image.*uold))/sum(sum(uold));
            c2=sum(sum(Image.*(1-uold)))/sum(sum(1-uold));

            I1=(Image-c1).^2;
            I2=(Image-c2).^2;
        end
        u_x=gradx(utilde);
        u_y=grady(utilde);
        
        % calcul de z_{k+1}
        z1=z1+sigma*u_x; z2=z2+sigma*u_y; 
        normeZ=norm_eps(z1,z2,0);
        nCond=normeZ>1;
        z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
    
        % calcul de u_{k+1}
        u=uold+tho*(div(z1,z2)-lambda*(I1-I2));
        u=min(max(u,0),1);
        utilde=u+theta*(u-uold);
        
        ub=u>mu;
        
        err(niter)=norm(u-uold);
        J1(niter)=compute_energy_smooth(ub,I1,I2,lambda,0); 
        J2(niter)=compute_energy(u,z1,z2,I1,I2,lambda);
    end
    err(isnan(err))=[];
    J1(isnan(J1))=[];
    J2(isnan(J1))=[];
end