%% Initialisation 
addpath(genpath('../'));
load('Images/Masks','Masks');

SelectImage={'Square.jpg','eight.tif','TumeurCerveaubis.png','TumeurCerveau.png','Poumon.png','CerveauDetail1.png'};
ImSelect=3;
SelectMethod={'Bande étroite','Histogrammes'};
MethSelect=2;

%% Valeur de lambda + nbins à tester + var booléennes
lambda=1.e1;
nbins=2;

Import=true;
Bruitage=false;
%% Paramètres numériques
itermax=200; 
mu=0.1; beta=0.5;
stop_1=1.e-8; stop_2=1.e-8;
visibility='on';

%% Lecture des images
Image=double(imread(SelectImage{ImSelect}));
Mask=Masks{ImSelect};
Background=Backgrounds{ImSelect};
%% Normalisation des images + ajout de bruit
Image=Image_Normalisation(Image,"2D");

if Bruitage==true
    bruit=0.05;
        
    Image=Image+bruit*randn(size(Image));
    Image=Image_Normalisation(Image,"2D");
    
end

% Importation des masques ou création de ces derniers
if ~Import
    Mask=roipoly(Image);
    Background=roipoly(Image);
    close;
end

if strcmp(SelectMethod{MethSelect},'Bande étroite')
    % Definition des paramètres numeriques 
    lambda=1.e-3; gamma=1;  %lambda : paramètre d'attache aux donnés  %gamma : paramètre de convexité 
    tho_u=0.5; tho_z=0.25;    %tho_u & tho_z : pas de descente pour les algo Forward-Backward
    beta=5; mu=0.1;           %beta : paramètre de la bande étroite     %mu : paramètre de seuilage 
    stop_1=0.005; stop_2=0.005;   %sigma : critère d'arrêt sur la décroissance de la fonctionnelle
    
    reset_band=10; % Nombre d'iterations de l'algo Primal-Dual avant de recalculer la level-set
    itermax=100;

    %% Algorithme de Chambol-Pock avec bande étroite
    tic
    [ub,J,niter]=bande_etroite_CP(Image,mask1,beta,tho_u,tho_z,lambda,mu,gamma,stop_1,stop_2,c1,c2,cinconnu,reset_band,itermax);
    t1=toc;
end

if strcmp(SelectMethod{MethSelect},'Histogrammes')
    fig=figure('Name','Histogrammes des régions à segmenter','NumberTitle','off');
    ax1=subplot(2,1,1); ax2=subplot(2,1,2);
    [g0,g1,T,sigma_1,sigma_2,sigma_3,b]=create_histo(Image,Mask,Background,nbins,ax1,ax2,visibility);
    
    tic
    [ub,J,err_u,err_J,niter]=histo_loco(double(Mask),g0,g1,b,T,sigma_1,sigma_2,sigma_3,lambda,beta,mu,stop_1,stop_2,itermax);
    toc
end

figure();
imagesc(Image); axis off; axis image;
colormap gray
hold on
contour(ub,'r','Linewidth',3);
title(['Resultat histogrammes pour \lambda= ',num2str(lambda),' et nbins= ',num2str(nbins)]);
hold off

%% Courbes de convergence
%semilogy(1:niter,CostFunc)
figure();
subplot(1,3,1)
plot(1:niter,J)
title('Fonction cout algo des histogrammes')
xlabel('iterations')
ylabel('J');
subplot(1,3,2)
plot(1:niter,err_u);
title('Condition de stagnation de la solution')
xlabel('iterations')
subplot(1,3,3)
plot(1:niter,err_u);
title('Condition de stagnation de la fonctionnelle')
xlabel('iterations')