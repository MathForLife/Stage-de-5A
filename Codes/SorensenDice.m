function SDC=SorensenDice(X,Y)
%% Fonction calculant l'indice de Sorensen-Dice, permettant de mesurer la ressemblance de 2 masques binaire
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
% X,Y : masques binaires
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT :
% SDC : Indice de Sorensen-Dice definit comme \frac{|X & Y|}{|X|+|Y|}

sx=size(X); sy=size(Y);
if ~isequal(sx,sy)
    fprintf('Dim %d : masque X = %d  masque Y = %d \n',[1:2;sx;sy]);
    error('Les masques ne sont pas de dimension egales');
else
    dim=length(sx);
    if dim==2
        SDC=sum(sum(X & Y,1),2)/sum(sum(X | Y,1),2);
    elseif dim==3
        SDC=sum(sum(sum(X & Y,1),2),3)/sum(sum(sum(X | Y,1),2),3);
    end 
end
end