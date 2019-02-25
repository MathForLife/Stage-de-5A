function M=gradx(I)
%Calcul le gradient en x d'une image I
%Syntaxe: gradx(I)

[m,n,p]=size(I);
M=zeros(m,n,p);

M(1:end-1,:,:)=I(2:end,:,:)-I(1:end-1,:,:);


