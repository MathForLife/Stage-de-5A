function [ul,J,niter]=histo_loco(mask,g0,g1,q2_0,q3_0,b,T,sigma_1,sigma_2,sigma_3,mu,beta,s,stop_1,stop_2,itermax)
    niter=1;
    
    % Initialisation de la bande �troite
    u=mask;
%     D=signed_distance_from_mask(mask);
%     Bl=abs(D)<=beta;
    %Bl=ones(size(Image));
    
    % Initialisation des variables u et u_tilde
%     u=ul.*Bl; u_tilde=ul.*Bl;
    % Initialisation de la variabe duale z
    z1=gradx(u); z2=grady(u);
    normeZ=norm_eps(z1,z2,0);
    nCond=normeZ>1;
    z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
    % Initialisation des 2 variables duales q2 et q3 associ�es � des
    % histogrammes 
    q2=q2_0; q3=q3_0;
    
    J(niter)=compute_energy_histo(u,z1,z2,g0,g1,mu,beta);
    %% Premi�re it�ration de l'algorithme
    niter=2; 
    
    u_old=u;
    
    u=u-T.*(Operator(g1,q2,true)-Operator(g0,q3,true)-div(z1,z2));
    u=min(max(u,0),1);
    
    u_tilde=2*u+u_old;
    u_x=gradx(u_tilde); u_y=grady(u_tilde);
    
    z1=(z1+sigma_1*u_x); z2=(z2+sigma_1*u_y); 
    % Projection sur la boule unit�
    normeZ=norm_eps(z1,z2,0);
    nCond=normeZ>1;
    z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
     
    q2=q2+diag(sigma_2)*Operator(g1,u_tilde,false);
    q2=min(max(q2,-mu/beta),mu/beta);
    
    q3=q3+diag(sigma_3)*(Operator(g0,u_tilde,false)-b);
    q3=min(max(q3,-mu/(1-beta)),mu/(1-beta));
    
    J(niter)=compute_energy_histo(u,z1,z2,g0,g1,mu,beta);
    
    %% Boucle while
    while (niter<itermax && abs(J(niter)-J(niter-1))>stop_1 && norm(u-u_old)>stop_2)
        niter=niter+1;
        
        u_old=u;
    
        u=u-T.*(Operator(g1,q2,true)-Operator(g0,q3,true)-div(z1,z2));
        u=min(max(u,0),1);

        u_tilde=2*u+u_old;
        u_x=gradx(u_tilde); u_y=grady(u_tilde);

        z1=(z1+sigma_1*u_x); z2=(z2+sigma_1*u_y); 
        % Projection sur la boule unit�
        normeZ=norm_eps(z1,z2,0);
        nCond=normeZ>1;
        z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);

        q2=q2+diag(sigma_2)*Operator(g1,u_tilde,false);
        q2=min(max(q2,-mu/beta),mu/beta);

        q3=q3+diag(sigma_3)*(Operator(g0,u_tilde,false)-b);
        q3=min(max(q3,-mu/(1-beta)),mu/(1-beta));

        J(niter)=compute_energy_histo(u,z1,z2,g0,g1,mu,beta);
%         if (mod(niter,nb)==0 || max(max(abs(u-ul)))>beta)
%             ul=u>mu;
%             %u=ul; u_tilde=ul; 
%             
%             D=signed_distance_from_mask(ul);
%             Bl=abs(D)<=beta;
%             %Bl=ones(size(Image));
%         end
    end
    ul=u>s;
    J(isnan(J))=[];
end