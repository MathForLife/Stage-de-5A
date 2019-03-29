function nabla_u=grad_mat(u,sz)
%% Calcul du gradient d'une image 2D ou 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
% u : masque défini sur \Omega, mais écrit sous forme d'un vecteur 
% Les dimension de l'image originale sont contenues dans le vecteur sz
% Si l'image est en 3D alors sz possède 3 coordonnées représantant les dimension de l'image selon chaque axe
% On suppose ne travailler qu'avec des images en niveau de gris
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% nabla_u : gradient du masque u calcule par differences finies avec conditions de Dirichlet 
% Chaque ligne est correspond a un point de l'espace et chaque colonne a une direction d'etude de la derivee
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim=length(sz);
vect_length=length(u);

u_temp=reshape(u,sz);
    
if dim==2
    ux=zeros(sz); uy=zeros(sz);
    
    ux(1:end-1,:)=u_temp(2:end,:)-u_temp(1:end-1,:);
    uy(:,1:end-1)=u_temp(:,2:end)-u_temp(:,1:end-1);
    
    nabla_u=[reshape(ux,vect_length,1),reshape(uy,vect_length,1)];
elseif dim==3
    ux=zeros(sz); uy=zeros(sz); uz=zeros(sz);
    
    ux(1:end-1,:,:)=u_temp(2:end,:,:)-u_temp(1:end-1,:,:);
    uy(:,1:end-1,:)=u_temp(:,2:end,:)-u_temp(:,1:end-1,:);
    uz(:,:,1:end-1)=u_temp(:,:,2:end)-u_temp(:,:,1:end-1);
    
    nabla_u=[reshape(ux,vect_length,1),reshape(uy,vect_length,1),reshape(uz,vect_length,1)];
end
end