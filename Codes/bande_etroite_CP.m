function [ul,J, niter]=bande_etroite_CP(Image,mask,beta,tho_u,tho_z,lambda,mu,gamma,stop_1,stop_2,c1,c2,cinconnu,nb,itermax)
    niter=1;
    
    % Initialiation de I1 et I2
    I1=(Image-c1).^2; I2=(Image-c2).^2;
    
    % Initialisation de la bande étroite
    ul=mask;
    D=signed_distance_from_mask(mask);
    Bl=abs(D)<=beta;
    %Bl=ones(size(Image));
    % Initialisation des variables u, u_tilde et z
    u=ul.*Bl; u_tilde=ul.*Bl;
    z1=sqrt(2)*ul; z2=-sqrt(2)*ul;
     
    J(niter)=compute_energy(u,z1,z2,I1,I2,lambda,ul,Bl,gamma);
    %% Première itération de l'algorithme
    niter=2; 
    
    u_tilde_old=u_tilde;
    u_x=gradx(u); u_y=grady(u);
    
    z1=(z1+tho_z*u_x).*Bl; z2=(z2+tho_z*u_y).*Bl; 
    % Projection sur la boule unité
    normeZ=norm_eps(z1,z2,0);
    nCond=normeZ>1;
    z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
    
    u_tilde=(u_tilde+tho_u*(div(z1,z2)-lambda*(I1-I2)+gamma*ul))/(1+gamma*tho_u);
    u_tilde=min(max(u_tilde,0),1).*Bl;
    
    theta=(1+2*gamma*tho_u)^(-1/2);
    tho_u=theta*tho_u;
    tho_z=tho_z/theta;
    
    u=(u_tilde+theta*(u_tilde-u_tilde_old)).*Bl;
    
    J(niter)=compute_energy(u,z1,z2,I1,I2,lambda,ul,Bl,gamma); 
    
    %% Boucle while
    while (niter<itermax && abs(J(niter)-J(niter-1))/J(niter)>stop_1 && norm(u_tilde-u_tilde_old)/norm(u_tilde)>stop_2)
        niter=niter+1;
        
        u_tilde_old=u_tilde;
        u_x=gradx(u); u_y=grady(u);

        if cinconnu
            c1=sum(sum(Image.*u))/sum(sum(u));
            c2=sum(sum(Image.*(1-u)))/sum(sum(1-u));

            I1=(Image-c1).^2;
            I2=(Image-c2).^2;
        end
    
        z1=(z1+tho_z*u_x).*Bl; z2=(z2+tho_z*u_y).*Bl; 
        % Projection sur la boule unité
        normeZ=norm_eps(z1,z2,0);
        nCond=normeZ>1;
        z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);

        u_tilde=(u_tilde+tho_u*(div(z1,z2)-lambda*(I1-I2)+gamma*ul))/(1+gamma*tho_u);
        u_tilde=min(max(u_tilde,0),1).*Bl;

        theta=(1+2*gamma*tho_u)^(-1/2);
        tho_u=theta*tho_u;
        tho_z=tho_z/theta;

        u=(u_tilde+theta*(u_tilde-u_tilde_old)).*Bl;

        J(niter)=compute_energy(u,z1,z2,I1,I2,lambda,ul,Bl,gamma); 
        
        if (mod(niter,nb)==0 || max(max(abs(u-ul)))>beta)
            ul=u>mu;
            %u=ul; u_tilde=ul; 
            
            D=signed_distance_from_mask(ul);
            Bl=abs(D)<=beta;
            %Bl=ones(size(Image));
        end
    end
    ul=u>mu;
    J(isnan(J))=[];
end