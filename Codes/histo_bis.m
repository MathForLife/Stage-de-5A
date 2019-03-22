function [ub,J, err_u,err_J,niter]=histo_bis(Image,Foreground,Background,Nbins,lambda,Parameters,PlotOptions,StopCondition)
%% Algorithme de segmentation par histogrammes inspiré des travaux de R.Yildizoglu, J-F Aujol, N.Papadakis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% Image: Image que l'on cherche à segmenter en 2 régions. L'image est en nuances de gris et peut être en 2D ou en 3D 
% Foreground: Masque représentant une partie de la zone que l'on souhaite
% segmenter. Il est possible de faire une union de masques
% Background: Masque représentant une partie du fond de l'image. Il est
% aussi possible de faire une union de masques
% Nbins: Nombre de bins que l'on souhaite prendre en compte dans les histogrammes
%     -1: nbins de l'histogramme h1
%     -2: nbins de l'histogramme h0
% lambda: Paramètre de controle entre attache aux données et régularisation sur la fonction J à minimiser
% Parameters: Ensemble de paramètres numériques 
%     -1 mu: Seuil utilisé pour le calcul du masque binaire u_b=H(u-\mu)
%     -2 beta: Approximation du ratio Taille de la zone à segmenter/Taille de l'image
%     -3 theta: Paramètre utilisé pour le calcul de \tilde{u} dans l'algo de Chambolle-Pock
%     -4 epsilon: Paramètre de régularisation des conditions sur Sigma_2, Sigma_3 et Tau (epsilon>0)
% PlotOptions: Contient des informations sur l'affichage des histogramme (passage par une fonction Matlab)
%     -1: axe de l'histogramme h1
%     -2: axe de l'histogramme h0
%     -3: visibilité de l'histogramme (on/off)
% StopCondition: Contient les conditions d'arrêt de l'algorithme
%     -1 itermax: Maximum d'itérations dans la boucle while
%     -2 stop_u: Condition sur la décroissance de la fonction de masquage u
%     -3 stop_J: Condition sur la décroissance de la fonctionnelle J
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
% ub: Masque binaire sur l'image de la région segmentée 
% err_u: Vecteur contenant l'évolution du premier critère d'arrêt au cours des itérations
% err_J: Vecteur contenant l'évolution du deuxième critère d'arrêt au cours des itérations
% niter: Nombre d'itérations total de l'algorithme

% Récupération des paramètres numériques
mu=Parameters(1); beta=Parameters(2);theta=Parameters(3); epsilon=Parameters(4);
% Création des histogrammes
[A,B,Tau,Sigma_1,Sigma_2,Sigma_3,b]=create_histo(Image,Foreground,Background,Nbins,epsilon,PlotOptions);

niter=1; itermax=StopCondition(1);
% Quantites de controle
J=nan(1,itermax);
err_u=nan(1,itermax);
err_J=nan(1,itermax);

    
sz=size(Image);
dim=length(sz);
if dim==2
    vect_length=sz(1)*sz(2);
elseif dim==3
    vect_length=sz(1)*sz(2)*sz(3);
end
% Initialisation de la variable u et vectorisation de l'espace \Omega
u=double(Foreground); u=reshape(u,vect_length,1);
% Initialisation de la variabe duale Q1=P_B(\nabla u)
Q1=grad_mat(u,sz);
normQ1=norm_grad(Q1); normCond=normQ1>1;
Q1(normCond)=Q1(normCond)./normQ1(normCond);

% Initialisation des variables duales Q2 et Q3 valant respectivement :
% Q2=P_[-\lambda/\beta,-\lambda/\beta](Au)
% Q3=P_[-\lambda/(1-\beta),-\lambda/(1-\beta)](B(1-u))

Q2=min(max(A*u,-lambda/beta),lambda/beta);
Q3=min(max(B*(1-u),-lambda/(1-beta)),lambda/(1-beta));

J(niter)=compute_energy_histo(u,A,B,lambda,beta,sz);

%% Premiere iteration de l'algorithme
niter=2;

u_old=u;
% Calcul sur la variable primale
u=u-Tau.*(-div_mat(Q1,sz)+A'*Q2-B'*Q3);
u=min(max(u,0),1);

u_tilde=u+theta*(u-u_old);
% Calcul de la première variable duale
Q1=Q1+Sigma_1*grad_mat(u_tilde,sz);
normQ1=norm_grad(Q1); normCond=normQ1>1;
Q1(normCond)=Q1(normCond)./normQ1(normCond);
% Calcul des 2 autres variables duales
Q2=Q2+Sigma_2.*(A*u_tilde);
Q2=min(max(Q2,-lambda/beta),lambda/beta);

Q3=Q3-Sigma_3.*(B*u_tilde-b);
Q3=min(max(Q3,-lambda/(1-beta)),lambda/(1-beta));

ub=u>mu;
J(niter)=compute_energy_histo(ub,A,B,lambda,beta,sz);
% Calcul des conditions d'arrêt (les normes utilisés sont des normes 2 sur
% tous les pixels de l'image
cond_u=norm(u-u_old)/norm(u_old);
cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));

err_u(1)=10; err_u(niter)=cond_u;
err_J(1)=10; err_J(niter)=cond_J;

%% Boucle while
while (niter<StopCondition(1) && cond_u>StopCondition(2) && cond_J>StopCondition(3))
    niter=niter+1;
    
    u_old=u;
    % Calcul sur la variable primale
    u=u-Tau.*(-div_mat(Q1,sz)+A'*Q2-B'*Q3);
    u=min(max(u,0),1);
    
    u_tilde=u+theta*(u-u_old);
    % Calcul de la première variable duale
    Q1=Q1+Sigma_1*grad_mat(u_tilde,sz);
    normQ1=norm_grad(Q1); normCond=normQ1>1;
    Q1(normCond)=Q1(normCond)./normQ1(normCond);
    % Calcul des 2 autres variables duales
    Q2=Q2+Sigma_2.*(A*u_tilde);
    Q2=min(max(Q2,-lambda/beta),lambda/beta);
    
    Q3=Q3-Sigma_3.*(B*u_tilde-b);
    Q3=min(max(Q3,-lambda/(1-beta)),lambda/(1-beta));
    
    ub=u>mu;
    J(niter)=compute_energy_histo(ub,A,B,lambda,beta,sz);
    % Calcul des conditions d'arrêt (les normes utilisés sont des normes 2 sur
    % tous les pixels de l'image
    cond_u=norm(u-u_old)/norm(u_old);
    cond_J=abs(J(niter)-J(niter-1))/abs(J(niter-1));
    
    err_u(niter)=cond_u; err_J(niter)=cond_J;
end
u=reshape(u,sz); ub=reshape(ub,sz);

J(isnan(J))=[];
err_u(isnan(err_u))=[];
err_J(isnan(err_J))=[];
end