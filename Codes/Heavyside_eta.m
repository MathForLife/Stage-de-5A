function Hn=Heavyside_eta(z,eta,n)
%% Fonction de Heaviside regularisee 
% INPUTS :
% z : vecteur (ou matrice) de discretisation de l'espace
% eta : parametre de lissage 
% n : entier indiquant le lissage a prendre en compte (peut prendre les valeurs 1 et 2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT :
% Hn : valeur de la fonction de Heaviside sur l'espace discretise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n==1
    Hn=zeros(size(z));
    
    cond=abs(z)<=eta;
    
    Hn(z>eta)=1;
    Hn(z<-eta)=0;
    Hn(cond)=0.5*(1+z(cond)/eta+sin(pi*z(cond)/eta)/pi);
elseif n==2
    Hn=0.5*(1+2/pi*atan(z/eta));
end
end