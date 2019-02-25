function J=compute_energy_histo(u,z1,z2,g0,g1,mu,beta,varargin)
    if nargin==7
        TV=sum(sum(u.*div(z1,z2)));
        J_h1=sum(abs(sum(sum(u.*g1,1),2)),3);
        J_h0=sum(abs(sum(sum((1-u).*g0,1),2)),3);
        J=TV+mu/beta*J_h1+mu/(1-beta)*J_h0;
    elseif nargin==10
        ul=varargin{1};
        Bl=varargin{2};
        gamma=varargin{3};
        
        TV=sum(sum(u.*div(z1,z2).*Bl));
        J_h1=sum(abs(sum(sum(u.*g1.*Bl,1),2)),3);
        J_h0=sum(abs(sum(sum((1-u).*g0.*Bl,1),2)),3);
        Dist_bande=sum(sum((u-ul).^2.*Bl));
        J=TV+mu/beta*J_h1+mu/(1-beta)*J_h0+gamma*Dist_bande;
    else 
        f = errordlg("Mauvais nombre d'arguments, choisir entre 6 ou 8",'Erreur en entrée');
        uniwait(f);
    end
   
end