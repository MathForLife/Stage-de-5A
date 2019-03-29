function n=norm_grad(nabla_u)
%% Calcul la norme d'un vecteur nabla_u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT :
% nabla_u : matrice dont les lignes correspondent aux points de l'espace etudie 
% et les colonnes correpondent aux valeur du champs dans les differentes directions de l'espace
% (ex : si nabla_u=grad_mat(u) alors la coordonnée (nabla_u)_{i,j} représente la j-ème dérivée 
% de la fonction u étudiée au point x_i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT : 
% n : vecteur colonne contenant la norme 2 du champ nabla_u en chaque point de l'espace etudie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim=size(nabla_u,2);
if dim==2
    n=sqrt(nabla_u(:,1).^2+nabla_u(:,2).^2);
elseif dim==3
    n=sqrt(nabla_u(:,1).^2+nabla_u(:,2).^2+nabla_u(:,3).^2);
end
end