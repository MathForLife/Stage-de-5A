 function [g0,g1,T,sigma_1,sigma_2,sigma_3,b]=create_histo(Image,mask1,mask0,nbins,ax1,ax2,visibility)
[m,n]=size(Image);

%subplot(2,nbIm,1)
h0=histogram(ax1,Image(mask0),'NumBins',nbins,'BinLimits',[0,1],'Normalization','probability','Visible',visibility);
title(ax1,'Histogramme du fond h^0');
%subplot(2,nbIm,2)
h1=histogram(ax2,Image(mask1),'NumBins',nbins,'BinLimits',[0,1],'Normalization','probability','Visible',visibility);
title(ax2,'Histogramme de la région à ségmenter : h^1');

g0=zeros(m,n,nbins); g1=zeros(m,n,nbins);
sigma_1=0.25; sigma_2=nan(nbins,1); sigma_3=nan(nbins,1);

for lambda= 1:nbins
    g0(:,:,lambda)= ((Image>=h0.BinEdges(lambda))&(Image<=h0.BinEdges(lambda+1)))-h0.Values(lambda);
    g1(:,:,lambda)= ((Image>=h1.BinEdges(lambda))&(Image<=h1.BinEdges(lambda+1)))-h1.Values(lambda);
    
    sigma_2(lambda)=1/(1+sum(sum(g1(:,:,lambda))));
    sigma_3(lambda)=1/(1+sum(sum(g0(:,:,lambda))));
end
T=1./(3+sum(abs(g0)+abs(g1),3)); 
b=sum(sum(g0,2),1); b=reshape(b,nbins,1);