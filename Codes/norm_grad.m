function n=norm_grad(nabla_u)
% Calcul la norme d'un vecteur nabla_u
%
% nabla_u est une matrice dont la i-ème colonnes est la coordonnées du champs de vecteur
% calculé sur tout l'espace
% ex : si nabla_u=grad_mat(u) alors la coordonnée (nabla_u)_{i,j} représente la j-ème dérivée 
% de la fonction u étudiée au point x_i. Dans ce cas, u a été vectorisée 
dim=size(nabla_u,2);
if dim==2
    n=sqrt(nabla_u(:,1).^2+nabla_u(:,2).^2);
elseif dim==3
    n=sqrt(nabla_u(:,1).^2+nabla_u(:,2).^2+nabla_u(:,3).^2);
end
end