function dn=delta_eta(z,eta,n)
%% Fonction de Dirac regularisee
% INPUTS :
% z : vecteur (ou matrice) de discretisation de l'espace
% eta : parametre de lissage
% n : entier indiquant le lissage a prendre en compte (peut prendre les valeurs 1 et 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT :
% dn : valeur de la fonction de Dirac sur l'espace discretise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n==1
    dn=zeros(size(z));
    
    cond=abs(z)<=eta;
    
    dn(cond)=0.5/eta*(1+cos(pi*z(cond)/eta));
elseif n==2
    dn=(eta*pi*(1+(z/eta).^2)).^-1;
end

end
