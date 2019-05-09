function J=compute_energy(u,I1,I0,lambda,epsilon,sz)
%% Calcul la fonction cout a minimiser
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% u: le masque binaire associe a l'image I
% Si I est de dimension [m,n] alors u est un vecteur de dimension [m*n,1]
% A: Operateur associe a la fonction g1 (matrice de dimension [Nbins(1),m*n]
% B: Operateur associe a la fonction g0 (matrice de dimension [Nbins(2),m*n]
% lambda: Paramètre permetant de contrôler l'attache aux donnees et la regularisation de la fonction
% mu: Estimation du ratio \frac{\Omega_1}{\Omega}
% sz: Liste contenant les dimensions de l'image de depart (et donc de u sous forme matricielle).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% J: Evaluation de la fonction cout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TV=sum(norm_grad(grad_mat(u,sz),epsilon));
Foreground=sum(I1*u);
Background=sum(I0*(1-u));
J=TV+0.5*lambda*Foreground+0.5*lambda*Background;
end