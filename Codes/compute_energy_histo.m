function J=compute_energy_histo(u,A,B,lambda,beta,sz)
%% Calcul la fonction coût à minimiser
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% u: le masque binaire associé à l'image I
% Si I est de dimension [m,n] alors u est un vecteur de dimension [m*n,1]
% A: Opérateur associé à la fonction g1 (matrice de dimension [Nbins(1),m*n]
% B: Opérateur associé à la fonction g0 (matrice de dimension [Nbins(2),m*n]
% lambda: Paramètre permetant de contrôler l'attache aux données et la régularisation de la fonction
% mu: Estimation du ratio \Omega_1/\Omega
% sz: Liste contenant les dimensions de l'image de départ (et donc de u
% sous forme matricielle).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% J: Evaluation de la fonction coût
%
TV=sum(norm_grad(grad_mat(u,sz)));
L1_A=sum(abs(A*u));
L1_B=sum(abs(B*(1-u)));
J=TV+lambda/beta*L1_A+lambda/(1-beta)*L1_B;
end