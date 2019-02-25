function [g0,g1,q2_0,q3_0,T,sigma_1,sigma_2,sigma_3,b]=create_histo(Image,mask1,mask2,nbins)
[m,n]=size(Image);

subplot(211)
h0=histogram(Image(mask1),nbins,'Normalization','pdf');
subplot(212)
h1=histogram(Image(mask2),nbins,'Normalization','pdf');

g0=zeros(m,n,nbins); g1=zeros(m,n,nbins);
sigma_1=0.2; sigma_2=nan(nbins,1); sigma_3=nan(nbins,1);
q2_0=nan(nbins,1); q3_0=nan(nbins,1);

for lambda= 1:nbins
    g0(:,:,lambda)= ((h0.BinEdges(lambda)<= Image)&(Image<h0.BinEdges(lambda+1)))-h0.Values(lambda);
    g1(:,:,lambda)= ((h1.BinEdges(lambda)<= Image)&(Image<h1.BinEdges(lambda+1)))-h1.Values(lambda);
    
    sigma_2(lambda)=1/(1+sum(sum(g1(:,:,lambda))));
    sigma_3(lambda)=1/(1+sum(sum(g0(:,:,lambda))));
    q2_0(lambda)=h0.Values(lambda);
    q3_0(lambda)=h1.Values(lambda);
end
T=1./(3+sum(abs(g0)+abs(g1),3)); 
b=sum(sum(g0,2),1); b=reshape(b,nbins,1);