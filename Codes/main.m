%% Initialisation 
addpath(genpath('../'));
SelectImage={'eight.tif','TumeurCerveaubis.png','Poumon.png','CerveauDetail1.png'};
Image=double(imread(SelectImage{1}));
% Image=zeros(248);
% Image(62:186,62:186)=255;

Selection={'Bande étroite','Histogrammes'};
%select='Bande étroite';
select='Histogrammes';

Image=Image_Normalisation(Image,"2D");

bruitage=false;
if bruitage==true
    bruit=0.1;
    Image=Image+bruit*255*randn(size(Image));
    
    Min=min(min(min(Image)),0);
    Max=max(max(max(Image)),255);
    
    Image=Image-Min;
    Image=Image*255/(Max-Min);
end
mask1=roipoly(Image);  c1=mean(Image(mask1)); 
mask2=roipoly(Image); close; c2=mean(Image(mask2));

if strcmp(select,Selection{1})
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

if strcmp(select,Selection{2})
    stop_1=0.001 ; stop_2=0.001;
    beta=0.5; nbins=10;
    mu=1.e2; s=0.9;
    itermax=100;
    
    [g0,g1,T,sigma_1,sigma_2,sigma_3,b]=create_histo(Image,mask1,mask2,nbins);
    
    tic
    [ub,J,niter]=histo_loco(mask1,g0,g1,b,T,sigma_1,sigma_2,sigma_3,mu,beta,s,stop_1,stop_2,itermax);
    t1=toc;
end

figure();
imagesc(Image); axis off; axis image;
colormap gray
hold on
contour(ub,'r','Linewidth',3);
title(['Resultat Chanbol et Pock pour ' num2str(niter),' iterations']);
sprintf('C1 = %5.1f, C2= %5.1f, niter= %d, t= %6.4f ',c1,c2,niter,t1)
hold off

%% Courbes de convergence
%semilogy(1:niter,CostFunc)
figure();
plot(1:niter,J)
title('Fonction cout algo Chambol-Pock')
xlabel('iterations')
