function M=gradz(I)
%Calcul le gradient en z d'une image I
%Syntaxe: grady(I)

[m,n,p]=size(I);
M=zeros(m,n,p);

M(:,:,1:end-1)=I(:,:,2:end)-I(:,:,1:end-1);
