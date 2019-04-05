function [Images, Foregrounds, Backgrounds,Regions,Gold_Standards]=ImportImageMasks(filename,extension,Im2Train,ImWithRegion,ChangeMasks,varargin)
%% Fonction permetant l'importation de masques et d'images existants ainsi que la modification des masques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation des dictionnaires
Images={}; Foregrounds={}; Backgrounds={}; Regions={}; Gold_Standards={};
%% Importation des images et des masques
fprintf('Importation des images et des masques\n');
for im=Im2Train
    Images{im}=imread([filename{im},extension]);
    load(['Foregrounds/',filename{im},'_FG.mat'],'FG');
    Foregrounds{im}=FG;
    load(['Backgrounds/',filename{im},'_BG.mat'],'BG');
    Backgrounds{im}=BG;
    if ismember(im,ImWithRegion)
        load(['Regions/',filename{im},'_Region.mat'],'R');
        Regions{im}=R;
    end    
    load(['Gold_Standards/',filename{im},'_GS.mat'],'GS');
    Gold_Standards{im}=GS;
end

%% Normalisation des images
fprintf('Normalisation des images\n');
for im=Im2Train
    Images{im}=Image_Normalisation(Images{im},"2D");
end

%% Creation de nouveaux masques 
if ChangeMasks
    fprintf('Creation de nouveaux masques\n');

    Foreground2Change=varargin{1}; Background2Change=varargin{2}; Region2Change=varargin{3};
    for f=Foreground2Change
        fprintf("Choisissez la region a segmenter sur l'image %s\n",filename{f});
        Foregrounds{f}=roipoly(Images{f});
    end
    close;
    
    for b=Background2Change
        fprintf("Choisissez un element de fond de l'image %d\n",filename{b});
        Backgrounds{b}=roipoly(Images{b});
    end
    close;
    
    for r=Region2Change
        fprintf("Choisissez la region de la radio a comparer avec le gold standard\n");
        Regions{r}=roipoly(Images{r});
    end
    close;
end
end