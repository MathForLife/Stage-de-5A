function SDC=SorensenDice(X,Y)
%% Fonction calculant l'indice de Sorensen-Dice, permettant de mesurer la ressemblance de 2 masques binaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
% X,Y : masques binaires
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT :
% SDC : Indice de Sorensen-Dice definit comme \frac{|X & Y|}{|X|+|Y|}

sx=size(X);
dim=length(sx);

if dim==2
    vect_length=sx(1)*sx(2);
elseif dim==3
    vect_length=sx(1)*sx(2)*sx(3);
end

X_temp=reshape(X,vect_length,1);
Y_temp=reshape(Y,vect_length,1);

SDC=2*sum(X_temp & Y_temp)/(sum(X_temp)+sum(Y_temp));
end