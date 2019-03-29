function n_eps=norm_eps(Ix,Iy,epsilon)
%% Calcul de la norme 2 pour une matrice avec un terme de lissage epsilon
n_eps=sqrt(Ix.^2+Iy.^2+epsilon.^2);
end