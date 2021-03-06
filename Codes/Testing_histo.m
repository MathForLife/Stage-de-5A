clear all; close all;

image_names={'BrainHole','BrainMeta_A','BrainTumor_A','BrainTumorDetail','Coins','chat',...
    'Flag','GBM_A','GeometricShape','Gliome003_S','Lung','Parenchyme_C','pneumopath_6_A','pneumopath_6_A_Detail','Square'};
extension='.png'; addpath(genpath('../Images/'));
%% Choix des images sur lesquelles entrainer les algos + importation et modification des masques
Im2Test=[4]; NbImages=length(Im2Test);
ImWithRegion=5:8;

import_masks=true; change_masks=false;  Bruitage=false; CompareTexture=true; 
import_textures=true; change_textures=false; choose_texture=true; verbose=true;
Foreground2Change=[11]; Background2Change=[11]; Region2Change=[]; Texture2Change=[];

[Images, Textures, Foregrounds, Backgrounds,~,~]=ImportImageMasks(image_names,extension,Im2Test,import_masks,ImWithRegion,change_masks,Foreground2Change,Background2Change,Region2Change,import_textures,choose_texture,change_textures,Texture2Change);

%% Valeur de lambda + nbins a tester
lambda=1.e0;
nbinsF=2; nbinsB=2;
Nbins=[nbinsF,nbinsB];
% Nombre de bins a prendre en compte pour la zone a segmenter (F=Front) et pour le fond de l'image (B=Back)


%% Parametres numeriques
itermax=200;
stop_u=1.e-8; stop_J=-1.e-8;
StopConditions=[itermax,stop_u,stop_J];

mu=0.9; beta=0.5; theta=1; epsilon=0.01; % parametre >0 servant a eviter des singularites dans la construction des histogrammes
Parameters=[mu, beta, theta, epsilon, verbose];

visibility='on'; cumulative.value=false;
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
f4=figure('Name','Images','NumberTitle','off');

%% Boucle for
iter=0;
for im=Im2Test
    
    iter=iter+1;
    [ub,J, err_u,err_J,nQ1,nQ2,nQ3,niter]=HistogrammeSegmentation(Images{im},lambda,Foregrounds{im},Backgrounds{im},Nbins,cumulative,visibility,Parameters,StopConditions,Textures{im});
    
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
    
    %% Evolution des quantités Q1, Q2 et Q3
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
    
    %% Affichage des masques initiaux
    figure(f4)
    subplot(1,NbImages,iter)
    imshow(Images{im});
    hold on
    contour(Foregrounds{im},'r','Linewidth',3);
    contour(Backgrounds{im},'g','Linewidth',3);
    hold off
    legend('\partial\Omega_1','\partial\Omega_0')
    title([image_names{im}, " and initial masks"]);
    
end

%% Sauvegarde des masques modifies
if change_masks
    SaveMasks(image_names, Foregrounds,Foreground2Change,Backgrounds, Background2Change,Regions,Region2Change);
end