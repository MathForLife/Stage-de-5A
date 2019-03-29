%% Initialisation 
clear all; close all;

addpath(genpath('../'));
load('Images/Masks','Masks');
load('Images/Backgrounds','Backgrounds');

SelectImage={'Square.jpg','eight.tif','TumeurCerveaubis.png','TumeurCerveau.png','Poumon.png','CerveauDetail3.png'};
ImSelect=3;
SelectMethod={'Bande etroite','Histogrammes'};
MethSelect=2;

%% Valeur de lambda + nbins a tester + var booleennes
lambda=1.e1;
nbinsF=4; nbinsB=4; 
% Nombre de bins a prendre en compte pour la zone a segmenter (F=Front) et pour le fond de l'image (B=Back)
Nbins=[nbinsF,nbinsB];

Import=true;
Bruitage=false;
%% Parametres numeriques
itermax=500; 
mu=0.5; beta=0.5; theta=0;
epsilon=1; % parametre >0 servant a eviter des singularites dans la construction des histogrammes
stop_u=-1.e-6; stop_J=-1.e-6;
visibility='on';

%% Lecture des images
Image=double(imread(SelectImage{ImSelect}));
Image=Image_Normalisation(Image,"2D");
% Importation des masques ou creation de ces derniers
if ~Import
    Mask=roipoly(Image);
    Background=roipoly(Image);
    close;
else
    Mask=Masks{ImSelect};
    Background=Backgrounds{ImSelect};
end
%% Normalisation des images + ajout de bruit
Image=Image_Normalisation(Image,"2D");

if Bruitage==true
    bruit=0.05;
        
    Image=Image+bruit*randn(size(Image));
    Image=Image_Normalisation(Image,"2D");
    
end

if strcmp(SelectMethod{MethSelect},'Bande etroite')
    % Definition des parametres numeriques 
    lambda=1.e-3; gamma=1;  %lambda : parametre d'attache aux donnes  %gamma : parametre de convexite 
    tho_u=0.5; tho_z=0.25;    %tho_u & tho_z : pas de descente pour les algo Forward-Backward
    beta=5; mu=0.1;           %beta : parametre de la bande etroite     %mu : parametre de seuilage 
    stop_u=0.005; stop_J=0.005;   %sigma : critere d'arret sur la decroissance de la fonctionnelle
    
    reset_band=10; % Nombre d'iterations de l'algo Primal-Dual avant de recalculer la level-set
    itermax=100;

    %% Algorithme de Chambol-Pock avec bande etroite
    tic
    [ub,J,niter]=bande_etroite_CP(Image,mask1,beta,tho_u,tho_z,lambda,mu,gamma,stop_u,stop_J,c1,c2,cinconnu,reset_band,itermax);
    t1=toc;
end

if strcmp(SelectMethod{MethSelect},'Histogrammes')
    fig=figure('Name','Histogrammes des regions e segmenter','NumberTitle','off');
    ax1=subplot(2,1,1); ax2=subplot(2,1,2);
    
    PlotOptions={ax1,ax2,visibility};
    StopConditions=[itermax,stop_u,stop_J];
    Parameters=[mu, beta, theta, epsilon];
    [ub,J, err_u,err_J,niter]=HistogrammeSegmentation(Image,Mask,Background,Nbins,lambda,Parameters,PlotOptions,StopConditions);
end

figure();
imagesc(Image); axis off; axis image;
colormap gray
hold on
contour(ub,'r','Linewidth',3);
title(['Resultat histogrammes pour \lambda= ',num2str(lambda),' et nbins= ',num2str(nbinsF)]);
hold off

%% Courbes de convergence
figure();
subplot(1,3,1)
plot(J)
title('Fonction cout algo des histogrammes')
xlabel('iterations')
ylabel('J');
subplot(1,3,2)
semilogy(err_u);
title('Condition de stagnation de la solution')
xlabel('iterations')
subplot(1,3,3)
plot(err_J);
title('Condition de stagnation de la fonctionnelle')
xlabel('iterations')
