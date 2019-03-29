function Energy=compute_energy_smooth(u,I1,I2,lambda,epsilon)
%% Calcul de la fonction cout avec un terme de lissage sur la variation totale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS: 
% u : masque binaire permetant la segmentation
% I1,I2 : erreur quadratique moyenne entre l'image et respectivement c1, c2
% lambda : parametre de regularisation de l'attache aux donnees
% epsilon : parametre de lissage du terme de variation totale
% OUTPUTS:
% Energy : valeur de la fonctionnelle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ux=gradx(u);
uy=grady(u);
n_eps=norm_eps(ux,uy,epsilon);

Energy=sum(sum(n_eps+lambda*I1.*u+lambda*I2.*(1-u)));
end