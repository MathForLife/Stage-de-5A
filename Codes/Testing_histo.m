clear all; close all;

image_names={'Square','GeometricShape','Coins','Flag','BrainTumor','BrainTumorDetail','BrainHole','Lung'};
texture_names={'Energy','Entropy','Correlation','IDM','Inertia','Cluster_Shade','Cluster_Prominence'};
extension='.png'; addpath(genpath('../Images/'));
%% Choix des images sur lesquelles entrainer les algos + importation et modification des masques
Im2Test=[6]; NbImages=length(Im2Test);
ImWithRegion=5:8;

ChangeMasks=false; Bruitage=false; CompareTexture=true; ImportTexture=false;
Foreground2Change=[1,2,3]; Background2Change=[1,2,3]; Region2Change=3; Texture2Change=[1,3,6];

[Images, Foregrounds, Backgrounds,~,~]=ImportImageMasks(image_names,extension,Im2Test,ImWithRegion,ChangeMasks,Foreground2Change,Background2Change,Region2Change);
if ImportTexture
    Textures=ImportTextures(image_names,Im2Select,texture_names,Text2Select);
else
    %Textures=struct([]);
    
    d_patch=5; d_glcm=1; ChooseTexture=true;
    for im=Im2Test
        
    %% Creation des cartes de texture
    Textures(im)=TextureMapping(Images{im},image_names{im},d_patch,d_glcm,ChooseTexture);
    end
end

%% Valeur de lambda + nbins a tester
lambda=1.e2;
nbinsF=5; nbinsB=5;
% Nombre de bins a prendre en compte pour la zone a segmenter (F=Front) et pour le fond de l'image (B=Back)
Nbins=[nbinsF,nbinsB];

%% Parametres numeriques
itermax=500;
stop_u=-1.e-8; stop_J=-1.e-8;
StopConditions=[itermax,stop_u,stop_J];

mu=0.5; beta=0.5; theta=1; epsilon=0.1; % parametre >0 servant a eviter des singularites dans la construction des histogrammes
Parameters=[mu, beta, theta, epsilon];

visibility=true; cumulative.value=false; 
if cumulative.value
    cumulative.normalisation='cdf';
else
    cumulative.normalisation='probability';
end

%% Ajout de bruit aux images
if Bruitage==true
    bruit=0.05;
    
    for im=Im2Test
        Images{im}=Images{im}+bruit*randn(size(Images{im}));
        Images{im}=Image_Normalisation(Images{im},"2D");
    end
end


%% Definition des figures
f1=figure('Name','Resultat de la segmentation par histogrammes','NumberTitle','off');
f2=figure('Name','Evolution des quantites de controle pour les histogrammes','NumberTitle','off');
f3=figure('Name','Evolution des variables duales','NumberTitle','off');

% fig=figure('Name','Histogrammes des regions a segmenter','NumberTitle','off');
% ax1=subplot(2,1,1); ax2=subplot(2,1,2);
% PlotOptions={ax1,ax2,visibility};

%% Boucle for
iter=0;
for im=Im2Test
    
    iter=iter+1;
    if CompareTexture
        [ub,J, err_u,err_J,nQ1,nQ2,nQ3,niter]=HistogrammeSegmentation(Images{im},lambda,Foregrounds{im},Backgrounds{im},Nbins,cumulative,visibility,CompareTexture,Parameters,StopConditions,Textures(im));
    else
        [ub,J, err_u,err_J,nQ1,nQ2,nQ3,niter]=HistogrammeSegmentation(Images{im},lambda,Foregrounds{im},Backgrounds{im},Nbins,cumulative,visibility,CompareTexture,Parameters,StopConditions);
    end
    figure(f1);
    subplot(1,NbImages,iter)
    imagesc(Images{im}); axis off; axis image;
    colormap gray
    hold on
    contour(ub,'r','Linewidth',3);
    title({image_names{im};['\lambda_H= ',num2str(lambda)]});
    hold off
        
    %% Courbes de convergence
    figure(f2);
    subplot(3,NbImages,iter)
    plot(J)
    title('Fonction cout algo des histogrammes');
    xlabel('iterations')
    ylabel('J');
    subplot(3,NbImages,iter+NbImages)
    semilogy(err_u);
    title('Condition de stagnation de la solution')
    xlabel('iterations')
    subplot(3,NbImages,iter+2*NbImages)
    plot(err_J);
    title('Condition de stagnation de la fonctionnelle')
    xlabel('iterations')
    
    %% Evolution des quantités Q2 et Q3
    figure(f3);
    subplot(3,NbImages,iter)
    plot(nQ1)
    title('norm(Q1)')
    xlabel('iterations')
    ylabel('Q1');
    subplot(3,NbImages,iter+NbImages)
    plot(nQ2)
    title('norm(Q2)')
    xlabel('iterations')
    ylabel('Q2');
    subplot(3,NbImages,iter+2*NbImages)
    plot(nQ3);
    title('norm(Q3)')
    xlabel('iterations')
    
end

%% Sauvegarde des masques modifies
if ChangeMasks
    SaveMasks(image_names, Foregrounds,Foreground2Change,Backgrounds, Background2Change,Regions,Region2Change);
end