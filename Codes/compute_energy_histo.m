function J=compute_energy_histo(u,A,B,lambda,beta,sz)
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
TV=sum(norm_grad(grad_mat(u,sz)));
L1_A=sum(abs(A*u));
L1_B=sum(abs(B*(1-u)));
J=TV+lambda/beta*L1_A+lambda/(1-beta)*L1_B;
end