function [Foregrounds, Backgrounds,Regions,Gold_Standards]=ImportMasks(image_names,extension,Im2Consider,ImWithRegion)
%% Initialisation des dictionnaires
Foregrounds={}; Backgrounds={}; Regions={}; Gold_Standards={};
%% Importation des masques 
fprintf('Importation des masques\n');
for im=Im2Consider
    %% Foregrounds
    FG=logical(imread(['Foregrounds/',image_names{im},'_FG',extension]));
    Foregrounds{im}=Image_Normalisation(FG,"2D");
    %% Backgrounds
    BG=logical(imread(['Backgrounds/',image_names{im},'_BG',extension]));
    Backgrounds{im}=Image_Normalisation(BG,"2D");
    %% Regions
    if ismember(im,ImWithRegion)
        R=logical(imread(['Regions/',image_names{im},'_Region',extension]));
        Regions{im}=Image_Normalisation(R,"2D");    
    end
    %% Gold Standars
    GS=logical(imread(['Gold_Standards/',image_names{im},'_GS',extension]));
    Gold_Standards{im}=Image_Normalisation(GS,"2D");
end