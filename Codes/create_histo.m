function [A,B,Tau,Sigma_1,Sigma_2,Sigma_3,b]=create_histo(Image,Foreground,Background,Nbins,Cumulative,epsilon,PlotOptions)
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
[m,n,p]=size(Image);

if Cumulative
    h1=histogram(PlotOptions{1},Image(Foreground),'NumBins',Nbins(1),'BinLimits',[0,1],'Normalization','cdf','Visible',PlotOptions{3});
    h0=histogram(PlotOptions{2},Image(Background),'NumBins',Nbins(2),'BinLimits',[0,1],'Normalization','cdf','Visible',PlotOptions{3});
else
    h1=histogram(PlotOptions{1},Image(Foreground),'NumBins',Nbins(1),'BinLimits',[0,1],'Normalization','probability','Visible',PlotOptions{3});
    h0=histogram(PlotOptions{2},Image(Background),'NumBins',Nbins(2),'BinLimits',[0,1],'Normalization','probability','Visible',PlotOptions{3});
end
title(PlotOptions{1},'Histogramme de la region a segmenter : h^1');
title(PlotOptions{2},'Histogramme du fond h^0');

%% Initialisation des operateurs A et B
if p>1
    A=zeros(Nbins(1),m,n,p); B=zeros(Nbins(2),m,n,p);
else
    A=zeros(Nbins(1),m,n); B=zeros(Nbins(2),m,n);
end

%% Boucle sur les bins de h1
for lambda= 1:Nbins(1)
    if Cumulative
        A(lambda,:,:,:)= double(Image<=h1.BinEdges(lambda+1))-h1.Values(lambda);
    else
        A(lambda,:,:,:)= double((Image>=h1.BinEdges(lambda))&(Image<=h1.BinEdges(lambda+1)))-h1.Values(lambda);
    end
    
end
%% Boucle sur les bins de h0
for lambda=1:Nbins(2)
    if Cumulative
        B(lambda,:,:,:)= double(Image<=h0.BinEdges(lambda+1))-h0.Values(lambda);
    else
        B(lambda,:,:,:)= double((Image>=h0.BinEdges(lambda))&(Image<=h0.BinEdges(lambda+1)))-h0.Values(lambda);
    end
end

%% Redimentionnement des matrices !
if p>1
    A=reshape(A,Nbins(1),m*n*p);
    B=reshape(B,Nbins(2),m*n*p);
else
    A=reshape(A,Nbins(1),m*n);
    B=reshape(B,Nbins(2),m*n);
end
%% Creation des matrices de regularisation
Sigma_1=1/(epsilon+4);
Sigma_2=1./(epsilon+sum(abs(A),2));
Sigma_3=1./(epsilon+sum(abs(B),2));

Tau=1./(2+sum(abs(A),1)+sum(abs(B),1)+epsilon)'; % On souhaite que Tau soit un vecteur de R^{m*n}
b=sum(B,2);