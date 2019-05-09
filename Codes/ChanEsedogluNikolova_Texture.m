function [ub, J, err_u,err_J, niter]=ChanEsedogluNikolova_Texture(Image, Textures,lambda, Foreground, Background, Cinconnu, Parameters, StopConditions)
%% Algorithme de Chan-Esedoglu-Nikolova calculant un masque binaire de segmentation
% INPUTS:
% Image: Image que l'on cherche à segmenter en 2 régions. L'image est en nuances de gris et peut etre en 2D ou en 3D 
% u0 : Masque représentant une partie de la zone que l'on souhaite segmenter. Il est possible de faire une union de masques
% lambda: Parametre de controle entre attache aux donnees et regularisation sur la fonction J à minimiser
% Parameters: Contient des informations propres a l'implementation numerique de CEN
%     -1 tho : Pas initial de l'algorithme de gradient projete
%     -2 mu : Seuil utilise pour le calcul du masque binaire u_b=H(u-\mu)
%     -3 epsilon :  Parametre de lissage pour le calcul de la norme 
% Colors: Contient des informations sur les parametres de couleurs donnes utilises pour la minimisation
%     -1 c1 : couleur moyenne de la region a segmenter
%     -2 c2 : couleur moyenne du fond de l'image
%     -3 cinconnu : variable booleenne indiquant si on modifie les couleurs
%     a chaque iteration (influe sur le caractere convexe de la fonction J)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debut algo
itermax=StopConditions(1); 
tho=Parameters(1); mu=Parameters(2); epsilon=Parameters(3);
%% Initialisation de toutes les variables
niter=1; 
% Variable interne permettant d'eviter les boucles infinies quand les conditions d'arret ne sont pas realisables
k=0; seuil=30; 
s=tho;

err_u=nan(1,itermax);
err_J=nan(1,itermax);
J=nan(1,itermax);

sz=size(Image); dim=length(sz);
if dim==3
    vect_length=sz(1)*sz(2)*sz(3);
else
    vect_length=sz(1)*sz(2);
end

[Indicator,c1,c0]=create_indicators(Textures, Foreground, Background, vect_length);

u=reshape(Image,vect_length,1); ub=u>mu;

I1=(Indicator-c1).^2; I0=(Indicator-c0).^2;

% Calcul des quantités de controle
J(1)=compute_energy(ub,I1,I0,lambda,epsilon,sz);
cond_u=10; cond_J=10; 
err_u(1)=cond_u; err_J(1)=cond_J;

%% Boucle while
while (niter<itermax && cond_u>StopConditions(2)&& cond_J>StopConditions(3))
    niter=niter+1;
    
    uold=u;
    
    if Cinconnu
        % Cas ou les couleurs sont aussi des variables de la fonctionnelle
        c1=Indicator*uold/sum(uold);
        c0=Indicator*(1-uold)/sum(1-uold);
        
        I1=(Indicator-c1).^2;
        I0=(Indicator-c0).^2;
    end
    % Calcul du gradient de u
    grad_u=grad_mat(u,sz);
    n_eps=norm_grad(grad_u,epsilon);
    
    % Descente de gradient projete
    u=uold+s*(div_mat(grad_u./n_eps,sz)-lambda*(sum(I1-I0,1)'));
    u=min(max(u,0),1);
    ub=u>mu;
    
    % Calcul des quantitees de controle
    J(niter)=compute_energy(ub,I1,I0,lambda,epsilon,sz);

    cond_u=norm(u-uold)/norm(uold);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    if J(niter)>J(niter-1)
        % Si la fonctionnelle augmente, on divise le pas par 2 et on recommence l'iteration
        fprintf('J= %f, niter= %d\n',J(niter),niter)
        niter=niter-1; k=k+1;
        
        u=uold;
        s=s/2;
    else
        % Sinon, on accepte la nouvelle iteration et on reinitialise le pas 
        s=tho; k=0;
        err_u(niter)=cond_u;
        err_J(niter)=cond_J;
    end
    
    if k>seuil
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