 function [g0,g1,Tau,Sigma_1,Sigma_2,Sigma_3,b]=create_histo(Image,Mask1,Mask0,Nbins,epsilon,PlotOptions)
%% Creation des quantités nécessaires à l'algorithme de segmentation par histogrammes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
% Image: Image en 2D ou en 3D et normalisée en nuances de gris sur laquelle on souhaite faire la segmentation 
% Mask1: Masque binaire d'une partie de la région que l'on souhaite segmenter
% Mask0: Masque binaire d'une partie du fond de l'image (peut être l'union de plusieurs masques
% Nbins: Nombre de bins que l'on souhaite prendre en compte dans les histogrammes
%     -1: nbins de l'histogramme h1
%     -2: nbins de l'histogramme h0
% epsilon: Paramètre de régularisation des conditions sur Sigma_2, Sigma_3 et Tau (epsilon>0)
% PlotOptions: Contient des informations sur l'affichage des histogramme (passage par une fonction Matlab)
%     -1: axe de l'histogramme h1
%     -2: axe de l'histogramme h0
%     -3: visibilité de l'histogramme (on/off)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS :
% g0, g1: quantités associées aux histogrammes h0 et h1 et définies comme g_i(x,\lambda)=Ind_{I(x)==\lambda}-h_i(\lambda)
% Les 2 opérateurs sont vectorisés en espace afin de faciliter les calculs. Ainsi, si l'image est de dimension [m,n]
%     -g1 est de dimension [nbins{1},m*n]
%     -g0 est de dimension [nbins{2},m*n]
% Tau: Vecteur de dimension [m*n,1] comportant une condition sur chaque coordonnée d'espace
% Sigma_i: Vecteur de dimension [nbins{i},1] représentant la matrice de conditionnement des variables duales
% b: Vecteur de dimension [nbins{2},1] trouvé en calculant le produit scalaire dans \Omega de g0 avec 1

[m,n]=size(Image);

h1=histogram(PlotOptions{1},Image(Mask1),'NumBins',Nbins(1),'BinLimits',[0,1],'Normalization','probability','Visible',PlotOptions{3});
title(PlotOptions{1},'Histogramme de la region a segmenter : h^1');

h0=histogram(PlotOptions{2},Image(Mask0),'NumBins',Nbins(2),'BinLimits',[0,1],'Normalization','probability','Visible',PlotOptions{3});
title(PlotOptions{2},'Histogramme du fond h^0');

g1=zeros(Nbins(1),m,n); g0=zeros(Nbins(2),m,n);
%sigma_1=0.25; sigma_2=nan(nbins,1); sigma_3=nan(nbins,1);

% Boucle sur les bins de h1
for lambda= 1:Nbins(1)
    g1(lambda,:,:)= double((Image>=h1.BinEdges(lambda))&(Image<=h1.BinEdges(lambda+1)))-h1.Values(lambda);
end
% Boucle sur les bins de h0
for lambda=1:Nbins(2)
    g0(lambda,:,:)= double((Image>=h0.BinEdges(lambda))&(Image<=h0.BinEdges(lambda+1)))-h0.Values(lambda);
    %     sigma_2(lambda)=1/(1+sum(sum(g1(:,:,lambda))));
%     sigma_3(lambda)=1/(1+sum(sum(g0(:,:,lambda))));
end
% T=1./(2+sum(abs(g0)+abs(g1),3)); 
% b=sum(sum(g0,2),1);

%% Redimentionnement des matrices !
g1=reshape(g1,Nbins(1),m*n);
g0=reshape(g0,Nbins(2),m*n);

Sigma_1=1/(epsilon+4);
Sigma_2=1./(epsilon+sum(g1,2));
Sigma_3=1./(epsilon+sum(g0,2));

Tau=1./(2+sum(abs(g1),1)+sum(abs(g0),1)+epsilon)'; % On souhaite que Tau soit un vecteur de R^{m*n} 
b=sum(g0,2);