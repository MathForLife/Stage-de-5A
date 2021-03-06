function [Images, Textures, Foregrounds, Backgrounds, Regions,Gold_Standards]=ImportImageMasks(image_names,extension,Im2Consider,import_masks,ImWithRegion,change_masks,Foreground2Change,Background2Change,Region2Change,import_textures,choose_texture,change_textures,Texture2Change)
% Fonction permetant l'importation de masques et d'images existants ainsi que la modification des masques
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Importation des images
Images=ImportImages(image_names,extension,Im2Consider); 
%% Importation ou creation des masques
Gold_Standards={};
if import_masks
    [Foregrounds, Backgrounds,Regions,Gold_Standards]=ImportMasks(image_names,extension,Im2Consider,ImWithRegion);
    %% Creation de nouveaux masques
    if change_masks
        [Foregrounds, Backgrounds,Regions]=CreatMasks(Images,image_names,extension,Foregrounds,Foreground2Change,Backgrounds,Background2Change,Regions,Region2Change);
    end
else
    [Foregrounds, Backgrounds,Regions]=CreatMasks(Images,image_names,extension,{},Im2Consider,{},Im2Consider,{},Im2Consider);
end

%% Importation ou creation des textures
if import_textures
    Textures=ImportTextures(Images,image_names,Im2Consider,choose_texture);
    %% Creation de nouvelles textures
    if change_textures
        Textures=CreatTextures(Images,Textures,image_names,Texture2Change,choose_texture);
    end
else
    Textures=CreatTextures(Images,{},image_names,Im2Consider,choose_texture);
end
    
end