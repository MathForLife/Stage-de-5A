function [ub,J, err_u,err_J,niter]=HistogrammeSegmentation(Image,Foreground,Background,Nbins,lambda,Parameters,PlotOptions,StopCondition)
%% Algorithme de segmentation par histogrammes inspire des travaux de R.Yildizoglu, J-F Aujol, N.Papadakis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% Image: Image que l'on cherche a segmenter en 2 regions. L'image est en nuances de gris et peut être en 2D ou en 3D 
% Foreground: Masque representant une partie de la zone que l'on souhaite segmenter. Il est possible de faire une union de masques
% Background: Masque representant une partie du fond de l'image. Il est aussi possible de faire une union de masques
% Nbins: Nombre de bins que l'on souhaite prendre en compte dans les histogrammes
%     -1: nbins de l'histogramme h1
%     -2: nbins de l'histogramme h0
% lambda: Paramètre de controle entre attache aux donnees et regularisation sur la fonction J a minimiser
% Parameters: Ensemble de paramètres numeriques 
%     -1 mu: Seuil utilise pour le calcul du masque binaire u_b=H(u-\mu)
%     -2 beta: Approximation du ratio Taille de la zone a segmenter/Taille de l'image
%     -3 theta: Paramètre utilise pour le calcul de \tilde{u} dans l'algo de Chambolle-Pock
%     -4 epsilon: Paramètre de regularisation des conditions sur Sigma_2, Sigma_3 et Tau (epsilon>0)
% PlotOptions: Contient des informations sur l'affichage des histogramme (passage par une fonction Matlab)
%     -1: axe de l'histogramme h1
%     -2: axe de l'histogramme h0
%     -3: visibilite de l'histogramme (on/off)
% StopCondition: Contient les conditions d'arrêt de l'algorithme
%     -1 itermax: Maximum d'iterations dans la boucle while
%     -2 stop_u: Condition sur la decroissance de la fonction de masquage u
%     -3 stop_J: Condition sur la decroissance de la fonctionnelle J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% ub: Masque binaire sur l'image de la region segmentee 
% err_u: Vecteur contenant l'evolution du premier critère d'arrêt au cours des iterations
% err_J: Vecteur contenant l'evolution du deuxième critère d'arrêt au cours des iterations
% niter: Nombre d'iterations total de l'algorithme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debut algo
% Recuperation des paramètres numeriques
itermax=StopCondition(1);
mu=Parameters(1); beta=Parameters(2);theta=Parameters(3); epsilon=Parameters(4);
% Creation des histogrammes
[A,B,Tau,Sigma_1,Sigma_2,Sigma_3,b]=create_histo(Image,Foreground,Background,Nbins,epsilon,PlotOptions);
% Calcul de la dimension de l'image pour sa vectorisation
sz=size(Image);
dim=length(sz);
if dim==2
    vect_length=sz(1)*sz(2);
elseif dim==3
    vect_length=sz(1)*sz(2)*sz(3);
end

%% Initialisation de toutes les variables
niter=1; 
% Variable interne permettant d'eviter les boucles infinies quand les conditions d'arret ne sont pas realisables
k=0; seuil=30; 
Tu=Tau; S1=Sigma_1; S2=Sigma_2; S3=Sigma_3;

J=nan(1,itermax);
err_u=nan(1,itermax);
err_J=nan(1,itermax);
    
% Initialisation de la variable u et vectorisation de l'espace \Omega
u=double(Foreground); u=reshape(u,vect_length,1);
v=double(Background); v=reshape(v,vect_length,1);

% Initialisation de la variabe duale Q1=P_B(\nabla u)
Q1=grad_mat(u,sz);
normQ1=norm_grad(Q1); normCond=normQ1>1;
Q1(normCond)=Q1(normCond)./normQ1(normCond);

% Initialisation des variables duales Q2 et Q3 valant respectivement :
% Q2=P_{[-\lambda/\beta,-\lambda/\beta]}(Au)
% Q3=P_{[-\lambda/(1-\beta),-\lambda/(1-\beta)]}(B(1-u))
Q2=min(max(A*u,-lambda/beta),lambda/beta);
Q3=min(max(B*v,-lambda/(1-beta)),lambda/(1-beta));

% Calcul des quantites de controle
J(1)=compute_energy_histo(u,A,B,lambda,beta,sz);
cond_u=10; err_u(1)=cond_u; 
cond_J=10; err_J(1)=cond_J;

%% Boucle while
while (niter<itermax && cond_u>StopCondition(2) && cond_J>StopCondition(3))
    niter=niter+1;
    
    u_old=u;
    % Iteration sur la variable primale
    figure(3); 
    subplot(141); imagesc(reshape(u,sz)); axis off; axis image;
    u=u+Tu.*(div_mat(Q1,sz)-A'*Q2+B'*Q3);
    subplot(142); imagesc(reshape(u,sz)); axis off; axis image;
    u=min(max(u,0),1);
    subplot(143); imagesc(reshape(u,sz)); axis off; axis image;
    ub=u>mu;
    
    u_tilde=u+theta*(u-u_old);
    subplot(144); imagesc(reshape(u_tilde,sz)); axis off; axis image;
    
    % Iteration sur la première variable duale
    Q1=Q1+S1*grad_mat(u_tilde,sz);
    normQ1=norm_grad(Q1); normCond=normQ1>1;
    Q1(normCond)=Q1(normCond)./normQ1(normCond);
    
    % Iteration sur les 2 autres variables duales
    Q2=Q2+S2.*(A*u_tilde);
    Q2=min(max(Q2,-lambda/beta),lambda/beta);
    
    Q3=Q3-S3.*(B*u_tilde-b);
    Q3=min(max(Q3,-lambda/(1-beta)),lambda/(1-beta));
    
    % Calcul des quantitees de controle (les normes utilises sont des normes 2 sur tous les pixels de l'image)
    J(niter)=compute_energy_histo(ub,A,B,lambda,beta,sz);
    cond_u=norm(u-u_old)/norm(u_old);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));

    if J(niter)>J(niter-1)
        % Si la fonctionnelle augmente, on divise le pas par 2 et on recommence l'iteration
        fprintf('J= %f, niter= %d\n',J(niter),niter)
        niter=niter-1; k=k+1;
        
        u=u_old;
        Tu=Tu/2; S1=S1*2; S2=S2*2; S3=S3*2;
    else
        % Sinon, on accepte la nouvelle iteration et on reinitialise les pas 
        Tu=Tau; S1=Sigma_1; S2=Sigma_2; S3=Sigma_3;
        k=0;
        
        err_u(niter)=cond_u;
        err_J(niter)=cond_J;
    end
    
    if k>seuil
        % Si on a refait plus de 'seuil' iterations sans faire decroitre la fonctionnelle, on s'arrete
        J(niter+1)=nan;
        niter=itermax;
    end
    
end

%% Redimentionnement du masque binair et des vecteurs de controle de l'algo
ub=reshape(ub,sz);

J(isnan(J))=[];
err_u(isnan(err_u))=[];
err_J(isnan(err_J))=[];
end