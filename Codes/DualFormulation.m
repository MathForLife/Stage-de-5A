function [ub, J, err_u,err_J, niter]=DualFormulation(Image, u0, lambda, mu, Parameters,Colors, StopConditions)
%% Algorithme de Chambolle-Pock calculant un masque binaire de segmentation a partir d'une image et d'un masque initial
% INPUTS:
% Image: Image que l'on cherche à segmenter en 2 régions. L'image est en nuances de gris et peut etre en 2D ou en 3D
% u0 : Masque représentant une partie de la zone que l'on souhaite segmenter. Il est possible de faire une union de masques
% lambda: Parametre de controle entre attache aux donnees et regularisation sur la fonction J à minimiser
% mu: Seuil utilise pour le calcul du masque binaire u_b=H(u-\mu)
% Parameters: Contient des informations propres a CP
%     -1 tho_u : Pas initial de la descente de gradient sur la variable primale
%     -2 tho_z : Pas initial de la montee de gradient sur la variable duale
%     -3 theta : Parametre utilise pour le calcul de \tilde{u}_k dans CP
% Colors: Contient des informations sur les parametres de couleurs donnes utilises pour la minimisation
%     -1 c1 : couleur moyenne de la region a segmenter
%     -2 c2 : couleur moyenne du fond de l'image
%     -3 cinconnu : variable booleenne indiquant si on modifie les couleurs a chaque iteration (influe sur le caractere convexe de la fonction J)
% StopConditions: Contient les conditions d'arrêt de l'algorithme
%     -1 itermax: Maximum d'itérations dans la boucle while
%     -2 stop_u: Condition sur la décroissance de la fonction de masquage u
%     -3 stop_J: Condition sur la décroissance de la fonctionnelle J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% ub: Masque binaire sur l'image de la région segmentée
% J : Vecteur contenant l'evolution de la fonction cout au cours des itérations
% err_u: Vecteur contenant l'evolution du premier critère d'arrêt au cours des itérations
% err_J: Vecteur contenant l'evolution du deuxième critère d'arrêt au cours des itérations
% niter: Nombre d'itérations total de l'algorithme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debut algo
itermax=StopConditions(1); c1=Colors(1); c2=Colors(2);
tho_u=Parameters(1); tho_z=Parameters(2); theta=Parameters(3);
%% Initialisation de toutes les variables
niter=1;
% Variable interne permettant d'eviter les boucles infinies quand les conditions d'arret ne sont pas realisables
k=0; seuil=50; 
u=u0; utilde=u0;

J=nan(1,itermax);
err_u=nan(1,itermax); err_J=nan(1,itermax);
su=tho_u; sz=tho_z; % Pas temporaires des algorithmes de gradient projetes 
% Initialisation de z
z1=gradx(utilde); z2=grady(utilde);
normeZ=norm_eps(z1,z2,0);
nCond=normeZ>1;

z1(nCond)=z1(nCond)./normeZ(nCond);
z2(nCond)=z2(nCond)./normeZ(nCond);
% Initialiation de I1 et I2
I1=(Image-c1).^2; I2=(Image-c2).^2;

% Calcul des quantites de controle
J(1)=compute_energy_smooth(u0,I1,I2,lambda,0);
cond_u=10; err_u(1)=cond_u; 
cond_J=10; err_J(1)=cond_J;
%% Boucle while
while (niter<itermax && cond_u>StopConditions(2)&& cond_J>StopConditions(3)) 
    niter=niter+1;
    
    uold=u; z1_old=z1; z2_old=z2;
    
    if Colors(3)
        % Cas ou les couleurs sont aussi des variables de la fonctionnelle
        c1=sum(sum(Image.*uold))/sum(sum(uold));
        c2=sum(sum(Image.*(1-uold)))/sum(sum(1-uold));
        
        I1=(Image-c1).^2;
        I2=(Image-c2).^2;
    end
    
    % Iteration sur la variable duale \textbf{z}
    z1=z1_old+sz*gradx(utilde); z2=z2_old+sz*grady(utilde);
    normeZ=norm_eps(z1,z2,0);
    nCond=normeZ>1;
    z1(nCond)=z1(nCond)./normeZ(nCond); z2(nCond)=z2(nCond)./normeZ(nCond);
    
    % Iteration sur la variable primale u + calcul de \tilde{u}
    u=uold+su*(div(z1,z2)-lambda*(I1-I2));
    u=min(max(u,0),1);
    utilde=u+theta*(u-uold);
    
    ub=u>mu; % Masque binaire utilise pour la segmentation
    
    % Calcul des quantitees de controle
    J(niter)=compute_energy_smooth(ub,I1,I2,lambda,0);
    
    cond_u=norm(u-uold)/norm(uold);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    if J(niter)>J(niter-1)
        % Si la fonctionnelle augmente, on divise le pas par 2 et on recommence l'iteration
        niter=niter-1; k=k+1;
        
        u=uold; z1=z1_old; z2=z2_old;
        su=su/2; sz=2*sz;  % ! A voir s'il ne faut pas diviser le pas de la variable duale par 2 ! 
    else
        % Sinon, on accepte la nouvelle iteration et on reinitialise le pas 
        su=tho_u; sz=tho_z; k=0;
        err_u(niter)=cond_u;
        err_J(niter)=cond_J;
    end
    
    if k>seuil
        % Si on a refait plus de 'seuil' iterations sans faire decroitre la fonctionnelle, on s'arrete
        J(niter+1)=nan;
        niter=itermax;
    end
end
%% Redimentionnement des vecteurs de controle de l'algo
err_u(isnan(err_u))=[];
err_J(isnan(err_J))=[];
J(isnan(J))=[];
end