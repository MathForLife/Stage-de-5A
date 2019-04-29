function [A,B,Tau,Sigma_1,Sigma_2,Sigma_3,b]=create_histo(Image,Foreground,Background,Nbins,epsilon,Cumulative,Visibility)
%% Creation des quantites necessaires a l'algorithme de segmentation par histogrammes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
% Image: Image en 2D ou en 3D et normalisee en nuances de gris sur laquelle on souhaite faire la segmentation
% Mask1: Masque binaire d'une partie de la region que l'on souhaite segmenter
% Mask0: Masque binaire d'une partie du fond de l'image (peut être l'union de plusieurs masques
% Nbins: Nombre de bins que l'on souhaite prendre en compte dans les histogrammes
%     -1: nbins de l'histogramme h1
%     -2: nbins de l'histogramme h0
% Cumulative : Booleen permettant de choisir la methode des histogrammes cumules ou la version normale
% epsilon: Paramètre de regularisation des conditions sur Sigma_2, Sigma_3 et Tau (epsilon>0)
% PlotOptions: Contient des informations sur l'affichage des histogramme (passage par une fonction Matlab)
%     -1: axe de l'histogramme h1
%     -2: axe de l'histogramme h0
%     -3: visibilite de l'histogramme (on/off)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS :
% g0, g1: quantites associees aux histogrammes h0 et h1 et definies comme g_i(x,\lambda)=Ind_{I(x)==\lambda}-h_i(\lambda)
% Les 2 operateurs sont vectorises en espace afin de faciliter les calculs. Ainsi, si l'image est de dimension [m,n]
%     -A est de dimension [nbins{1},m*n]
%     -B est de dimension [nbins{2},m*n]
% Tau: Vecteur de dimension [m*n,1] comportant une condition sur chaque coordonnee d'espace
% Sigma_i: Vecteur de dimension [nbins{i},1] representant la matrice de conditionnement des variables duales
% b: Vecteur de dimension [nbins{2},1] trouve en calculant le produit scalaire dans \Omega de g0 avec 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debut algo
fprintf('Creation des histogrammes et des operateurs associes');
[m,n,p]=size(Image); 

%% Initialisation des parametres d'affichage
fig=figure('Name','Histogrammes des regions a segmenter','NumberTitle','off');

ax1=subplot(2,1,1); title('Histogramme de la region a segmenter : h^1');
ax2=subplot(2,1,2); title("Histogramme du fond de l'image : h^0");

%% Creation des histogrammes de couleur
h1=histogram(ax1,Image(Foreground),'NumBins',Nbins(1),'BinLimits',[0,1],'Normalization',Cumulative.normalisation,'Visible',Visibility);
h0=histogram(ax2,Image(Background),'NumBins',Nbins(2),'BinLimits',[0,1],'Normalization',Cumulative.normalisation,'Visible',Visibility);

%% Initialisation des operateurs A et B
A=zeros(Nbins(1),m,n,p); B=zeros(Nbins(2),m,n,p);

%% Boucle sur les bins de h1
% figure(20)
for lambda= 1:Nbins(1)
    if Cumulative.value
        A(lambda,:,:,:)= double(Image<=h1.BinEdges(lambda+1))-h1.Values(lambda);
    else
        A(lambda,:,:,:)= double((Image>=h1.BinEdges(lambda))&(Image<=h1.BinEdges(lambda+1)))-h1.Values(lambda);
    end
%     subplot(1,Nbins(1),lambda)
%     imagesc(reshape(A(lambda,:,:),m,n,p)); axis off; axis image;
end
%% Boucle sur les bins de h0
% figure(30)
for lambda=1:Nbins(2)
    if Cumulative.value
        B(lambda,:,:,:)= double(Image<=h0.BinEdges(lambda+1))-h0.Values(lambda);
    else
        B(lambda,:,:,:)= double((Image>=h0.BinEdges(lambda))&(Image<=h0.BinEdges(lambda+1)))-h0.Values(lambda);
    end    
%     subplot(1,Nbins(2),lambda)
%     imagesc(reshape(B(lambda,:,:,:),m,n,p)); axis off; axis image;
end

%% Redimentionnement des matrices !
A=reshape(A,Nbins(1),m*n*p);
B=reshape(B,Nbins(2),m*n*p);

%% Creation des matrices de regularisation
Sigma_1=1/(epsilon+4);
Sigma_2=1./(epsilon+sum(abs(A),2));
Sigma_3=1./(epsilon+sum(abs(B),2));

Tau=1./(2+sum(abs(A),1)+sum(abs(B),1)+epsilon)'; % On souhaite que Tau soit un vecteur de R^{m*n}
b=sum(B,2);

if strcmp(Visibility,'off') 
    close(fig);
end