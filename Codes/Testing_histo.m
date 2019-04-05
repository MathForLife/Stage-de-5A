clear all; close all;

filename={'Square','GeometricShape','Coins','BrainTumor','BrainTumorDetail','BrainHole','Lung'};
extension='.png'; addpath(genpath('../Images/'));
%% Choix des images sur lesquelles entrainer les algos + importation et modification des masques
Im2Test=[5]; NbImages=length(Im2Test);
ImWithRegion=4:7;

ChangeMasks=false; Bruitage=false; 
Foreground2Change=[1,2,3]; Background2Change=[1,2,3]; Region2Change=3;

[Images, Foregrounds, Backgrounds,~,~]=ImportImageMasks(filename,extension,Im2Test,ImWithRegion,ChangeMasks,Foreground2Change,Background2Change,Region2Change);

%% Valeur de lambda + nbins a tester
lambda=1.e2;
nbinsF=10; nbinsB=10;
% Nombre de bins a prendre en compte pour la zone a segmenter (F=Front) et pour le fond de l'image (B=Back)
Nbins=[nbinsF,nbinsB];

%% Parametres numeriques
itermax=500;
stop_u=-1.e-8; stop_J=-1.e-8;
StopConditions=[itermax,stop_u,stop_J];

mu=0.5; beta=0.5; theta=0; epsilon=1; % parametre >0 servant a eviter des singularites dans la construction des histogrammes
Parameters=[mu, beta, theta, epsilon];

visibility='on'; cumulative=true;

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

fig=figure('Name','Histogrammes des regions a segmenter','NumberTitle','off');
ax1=subplot(2,1,1); ax2=subplot(2,1,2);
PlotOptions={ax1,ax2,visibility};

%% Boucle for
iter=0;
for im=Im2Test
    iter=iter+1;
    [ub,J, err_u,err_J,niter]=HistogrammeSegmentation(Images{im},Foregrounds{im},Backgrounds{im},Nbins,cumulative,lambda,Parameters,PlotOptions,StopConditions);
    
    figure(f1);
    subplot(1,NbImages,iter)
    imagesc(Images{im}); axis off; axis image;
    colormap gray
    hold on
    contour(ub,'r','Linewidth',3);
    title({filename{im};['\lambda_H= ',num2str(lambda)]});
    hold off
        
    %% Courbes de convergence
    figure(f2);
    subplot(3,NbImages,iter)
    plot(J)
    title('Fonction cout algo des histogrammes')
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
end

%% Sauvegarde des masques modifies
if ChangeMasks
    SaveMasks(filename, Foregrounds,Foreground2Change,Backgrounds, Background2Change,Regions,Region2Change);
end