function J=compute_energy(u,z1,z2,I1,I2,lambda,varargin)
    if nargin==6
        J=sum(sum(u.*div(z1,z2)+lambda*(I1.*u+I2.*(1-u))));
    elseif nargin==9
        ul=varargin{1};
        Bl=varargin{2};
        gamma=varargin{3};
        
        J=sum(sum((u.*div(z1,z2)+lambda*(I1.*u+I2.*(1-u))+0.5*gamma*(u-ul).^2)));
    else 
        f = errordlg("Mauvais nombre d'arguments, choisir entre 6 ou 8",'Erreur en entrée');
        uniwait(f);
    end
        
end