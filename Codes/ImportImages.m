function Images=ImportImages(image_names,extension,Im2Consider)
%% Initialisation du dictionnaire
Images={};
%% Importation des images
fprintf('Importation des images\n');
for im=Im2Consider
    Image=double(imread([image_names{im},extension]));
    % Normalisation des images
    Image=Image_Normalisation(Image,"2D"); 
    Images{im}=Image;
end
