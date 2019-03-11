function [ub,J, err_u,err_J,niter]=histo_loco(mask,g0,g1,b,T,sigma_1,sigma_2,sigma_3,mu,s,beta,stop_1,stop_2,itermax)
    niter=1;
    % Quantit�s de contr�le
    J=nan(1,itermax);
    err_u=nan(1,itermax);
    err_J=nan(1,itermax);
    
    % Initialisation des variables u et u_tilde
    u=mask;
    % Initialisation de la variabe duale z
    z1=gradx(u); z2=grady(u); 
    normeZ=norm_eps(z1,z2,0);
    nCond=normeZ>1;
    z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
    
    % Initialisation des 2 variables duales q2 et q3 en projetant les
    % r�sultats Au_0 et B(1-u_0) sur les domaine associ�s � q2 et q3
    q2=min(max(Operator(g1,u,false),-mu/beta),mu/beta);
    q3=min(max(Operator(g0,1-u,false),-mu/(1-beta)),mu/(1-beta));
    
    J(niter)=compute_energy_histo(u,g0,g1,mu,beta);
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
    
    ub=u>s;
    J(niter)=compute_energy_histo(ub,g0,g1,mu,beta);
    % Calcul des erreurs
    cond_u=norm(u-u_old)/norm(u_old);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    err_u(1)=10; err_u(niter)=cond_u;
    err_J(1)=10; err_J(niter)=cond_J;
    
    %% Boucle while
    while (niter<itermax && cond_u>stop_1 && cond_J>stop_2)
        niter=niter+1;
        
        u_old=u;
    
        u=u-T.*(Operator(g1,q2,true)-Operator(g0,q3,true)-div(z1,z2));
        u=min(max(u,0),1);

        u_tilde=2*u-u_old;
        u_x=gradx(u_tilde); u_y=grady(u_tilde);

        z1=z1+sigma_1*u_x; z2=z2+sigma_1*u_y; 
        % Projection sur la boule unit�
        normeZ=norm_eps(z1,z2,0);
        nCond=normeZ>1;
        z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);

        q2=q2+diag(sigma_2)*Operator(g1,u_tilde,false);
        q2=min(max(q2,-mu/beta),mu/beta);

        q3=q3+diag(sigma_3)*(Operator(g0,u_tilde,false)-b);
        q3=min(max(q3,-mu/(1-beta)),mu/(1-beta));
        
        ub=u>s;
        J(niter)=compute_energy_histo(ub,g0,g1,mu,beta);
        % Calcul des erreurs
        cond_u=norm(u-u_old)/norm(u_old);
        cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));

        err_u(niter)=cond_u;
        err_J(niter)=cond_J;
    end
    J(isnan(J))=[];
    err_u(isnan(err_u))=[];
    err_J(isnan(err_J))=[];
end