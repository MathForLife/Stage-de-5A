function [ub, J, err_u,err_J, niter]=DualFormulation_Texture(Image, Textures, lambda, Foreground, Background, Cinconnu, Parameters, StopConditions)
%% Algorithme de Chambolle-Pock calculant un masque binaire de segmentation a partir d'une image et d'un masque initial
% INPUTS:
% Image: Image que l'on cherche à segmenter en 2 régions. L'image est en nuances de gris et peut etre en 2D ou en 3D
% u0 : Masque représentant une partie de la zone que l'on souhaite segmenter. Il est possible de faire une union de masques
% lambda: Parametre de controle entre attache aux donnees et regularisation sur la fonction J à minimiser
% Parameters: Contient des informations propres a CP
%     -1 tho_u : Pas initial de la descente de gradient sur la variable primale
%     -2 tho_z : Pas initial de la montee de gradient sur la variable duale
%     -3 mu : Seuil utilise pour le calcul du masque binaire u_b=H(u-\mu)
%     -4 theta : Parametre utilise pour le calcul de \tilde{u}_k dans CP
% Colors: Contient des informations sur les parametres de couleurs donnes utilises pour la minimisation
%     -1 c1 : couleur moyenne de la region a segmenter
%     -2 c2 : couleur moyenne du fond de l'image
%     -3 cinconnu : variable booleenne indiquant si on modifie les couleurs a chaque iteration (influe sur le caractere convexe de la fonction J)
% StopConditions: Contient les conditions d'arrêt de l'algorithme
%     -1 itermax: Maximum d'itérations dans la boucle while
%     -2 stop_u: Condition sur la décroissance de la fonction de masquage u
%     -3 stop_J: Condition sur la décroissance de la fonctionnelle J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% ub: Masque binaire sur l'image de la région segmentée
% J : Vecteur contenant l'evolution de la fonction cout au cours des itérations
% err_u: Vecteur contenant l'evolution du premier critère d'arrêt au cours des itérations
% err_J: Vecteur contenant l'evolution du deuxième critère d'arrêt au cours des itérations
% niter: Nombre d'itérations total de l'algorithme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debut algo
itermax=StopConditions(1); 
tho_u=Parameters(1); tho_z=Parameters(2); mu=Parameters(3); theta=Parameters(4);
%% Initialisation de toutes les variables
niter=1;
% Variable interne permettant d'eviter les boucles infinies quand les conditions d'arret ne sont pas realisables
k=0; seuil=50; 

J=nan(1,itermax);
err_u=nan(1,itermax); err_J=nan(1,itermax);
tu=tho_u; tz=tho_z; % Pas temporaires des algorithmes de gradient projetes 


sz=size(Image); dim=length(sz);
if dim==3
    vect_length=sz(1)*sz(2)*sz(3);
else
    vect_length=sz(1)*sz(2);
end

[Indicator,c1,c0]=create_indicators(Textures, Foreground, Background, vect_length);

u=reshape(Image,vect_length,1); utilde=u;

% Initialisation de z
Z=grad_mat(utilde,sz);
normeZ=norm_grad(Z,0);
nCond=normeZ>1;

Z(nCond,:)=Z(nCond,:)./normeZ(nCond);
% Initialiation de I1 et I2
I1=(Indicator-c1).^2; I0=(Indicator-c0).^2;

% Calcul des quantites de controle
J(1)=compute_energy(u,I1,I0,lambda,0,sz);
cond_u=10; err_u(1)=cond_u; 
cond_J=10; err_J(1)=cond_J;
%% Boucle while
while (niter<itermax && cond_u>StopConditions(2)&& cond_J>StopConditions(3)) 
    niter=niter+1;
    
    uold=u; Z_old=Z; 
    
    if Cinconnu
        % Cas ou les couleurs sont cree aussi des variables de la fonctionnelle
        c1=Indicator*uold/sum(uold);
        c0=Indicator*(1-uold)/sum(1-uold);
        
        I1=(Indicator-c1).^2;
        I0=(Indicator-c0).^2;
    end
    
    % Iteration sur la variable primale u + calcul de \tilde{u}
    u=uold+tu*(div_mat(Z,sz)-lambda*(sum(I1-I0,1)'));
    u=min(max(u,0),1);
    utilde=u+theta*(u-uold);
    
    ub=u>mu; % Masque binaire utilise pour la segmentation
    
    % Iteration sur la variable duale \textbf{z}
    Z=Z_old+tz*grad_mat(utilde,sz); 
    normeZ=norm_grad(Z,0);
    nCond=normeZ>1;
    
    Z(nCond,:)=Z(nCond,:)./normeZ(nCond); 
    % Calcul des quantitees de controle
    J(niter)=compute_energy(ub,I1,I0,lambda,0,sz);
    
    cond_u=norm(u-uold)/norm(uold);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    if J(niter)>2*J(niter-1)
        % Si la fonctionnelle augmente, on divise le pas par 2 et on recommence l'iteration
        fprintf('J= %f, niter= %d\n',J(niter),niter)
        niter=niter-1; k=k+1;
        
        u=uold; Z=Z_old; 
        tu=tu/2; tz=2*tz;  % ! A voir s'il ne faut pas diviser le pas de la variable duale par 2 ! 
    else
        % Sinon, on accepte la nouvelle iteration et on reinitialise le pas 
        tu=tho_u; tz=tho_z; k=0;
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

ub=reshape(ub,sz);
end