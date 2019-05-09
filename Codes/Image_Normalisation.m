function Image=Image_Normalisation(Image,Dimension)
%% Normalisation des niveaux de gris d'une image 2D ou 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
% Image : Image a normaliser. Ca peut etre une image 2D en couleur ou en niveau de gris ainsi qu'une image 3D
% Dimension : chaine de caractere indiquant si l'image est en 2D ou en 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT :
% Image : Image 2D ou 3D normalisee en niveau de gris dont chaque pixel est a valeur dans [0,1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,p]=size(Image);

if strcmp(Dimension,"2D") && p>1
    if isa(Image,'logical')
        Image=Image(:,:,1);
    elseif isa(Image,'double')
        Image=mean(Image,3);
    end
end

if isa(Image,'double')
    Min=min(min(Image));
    Max=max(max(Image));

    Image=Image-Min;
    Image=Image/(Max-Min);
end

end