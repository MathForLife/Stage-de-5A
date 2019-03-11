function J=compute_energy_histo(u,g0,g1,mu,beta,varargin)
    if nargin==5
        TV=sum(sum(norm_eps(gradx(u),grady(u),0)));
        J_h1=sum(abs(sum(sum(u.*g1,1),2)),3);
        J_h0=sum(abs(sum(sum((1-u).*g0,1),2)),3);
        J=TV+mu/beta*J_h1+mu/(1-beta)*J_h0;
    else 
        f = errordlg("Mauvais nombre d'arguments, choisir entre 6 ou 8",'Erreur en entrée');
        uniwait(f);
    end
   
end