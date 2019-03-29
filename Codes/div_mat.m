function dvu=div_mat(nabla_u,sz)
%% Calcul la divergence d'un vecteur sur l'ensemble de l'espace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS : 
% nabla_u : matrice dont la i-ème colonnes est la coordonnées du champs de vecteur
% calculé sur tout l'espace
% ex : si nabla_u=grad_mat(u) alors la coordonnée (nabla_u)_{i,j} représente la j-ème dérivée 
% de la fonction u étudiée au point x_i. Dans ce cas, u a été vectorisée 
% 
% sz : vecteur contenant les dimensions de la fonction u avant qu'elle n'ai été vectorisée
% (ex : si u est définie sur une image 3D alors sz=[nx,ny,nz])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% dvu : vecteur de meme dimension que u representant la divergence du champ \nabla u
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Debut de l'algo
dim=length(sz);
vect_length=length(nabla_u(:,1)); % Calcul de nx*ny*nz

dvu_temp=zeros(sz);

%% Cas des images 2D
if dim==2
    ux_temp=reshape(nabla_u(:,1),sz); uy_temp=reshape(nabla_u(:,2),sz);

    Mx=dvu_temp; My=dvu_temp;

    Mx(2:end-1,:)=ux_temp(2:end-1,:)-ux_temp(1:end-2,:);
    Mx(1,:)=ux_temp(1,:);
    Mx(end,:)=-ux_temp(end-1,:);

    My(:,2:end-1)=uy_temp(:,2:end-1)-uy_temp(:,1:end-2);
    My(:,1)=uy_temp(:,1);
    My(:,end)=-uy_temp(:,end-1);

    dvu_temp=Mx+My;
%% Cas des images 3D
elseif dim==3
    ux_temp=reshape(nabla_u(:,1),sz); uy_temp=reshape(nabla_u(:,2),sz); uz_temp=reshape(nabla_u(:,3),sz);
    
    Mx=dvu_temp; My=dvu_temp; Mz=dvu_temp;

    Mx(2:end-1,:,:)=ux_temp(2:end-1,:,:)-ux_temp(1:end-2,:,:);
    Mx(1,:,:)=ux_temp(1,:,:);
    Mx(end,:,:)=-ux_temp(end-1,:,:);

    My(:,2:end-1,:)=uy_temp(:,2:end-1,:)-uy_temp(:,1:end-2,:);
    My(:,1,:)=uy_temp(:,1,:);
    My(:,end,:)=-uy_temp(:,end-1,:);

    Mz(:,:,2:end-1)=uz_temp(:,:,2:end-1)-uz_temp(:,:,1:end-2);
    Mz(:,:,1)=uz_temp(:,:,1);
    Mz(:,:,end)=-uz_temp(:,:,end-1);

    dvu_temp=Mx+My+Mz;
end
%% Vectorisation du resultat
dvu=reshape(dvu_temp,vect_length,1);
end