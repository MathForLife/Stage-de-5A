addpath(genpath('../'))
load('Images/Masks','Masks');
load('Images/Backgrounds','Backgrounds');

%% Valeur de lambda a tester + var booleennes
lambda_CP=100; lambda_CEN=100; %tester 6 et 1 pour 
Im2Test=[1,2,3];

Import=true;
Bruitage=false;

%% Parametres numeriques
itermax=200; cinconnu=false;
mu=0.1; eps=1;
theta=1;
tho_u=0.5; tho_z=0.25;
stop_1=1.e-8; stop_2=1.e-8;
hist_visibility='on';

%% Lecture des images
Image1=double(imread('Square.jpg'));
Image2=double(imread('eight.tif'));
Image3=double(imread('TumeurCerveaubis.png'));

NbImages=length(Im2Test);
Images={Image1,Image2,Image3};
%% Normalisation des images
for i=Im2Test
    Images{i}=Image_Normalisation(Images{i},"2D");
end
%% Ajout de bruit aux images
if Bruitage==true
    bruit=0.05;
    
    for i=Im2Test
        Images{i}=Images{i}+bruit*randn(size(Images{i}));
        Images{i}=Image_Normalisation(Images{i},"2D");
    end
end

% Importation des masques ou creation de ces derniers
if ~Import
    for m=Mask2Change
        Masks{m}=roipoly(Images{m});
    end
    close;
    for b=Background2Change
        Backgrounds{b}=roipoly(Images{b});
    end
    close;
end

%% Definition des figures
f1=figure('Name','Resultat de la segmentation','NumberTitle','off');
f2=figure('Name','Evolution des quantites de controle pour C-P','NumberTitle','off');
f3=figure('Name','Evolution des quantites de controle pour C-E-N','NumberTitle','off');

%% Boucle sur les images
iter=0;
for i=Im2Test
    iter=iter+1;
    
    c1=mean(Images{i}(Masks{i})); c2=mean(Images{i}(~Masks{i}));
        
    [Ub_CP, J_CP, Err_u_Cp,Err_J_CP, Niter_CP]=DualFormulation(Images{i}, double(Masks{i}), lambda_CP, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax);
    [Ub_CEN, J_CEN, Err_u_CEN,Err_J_CEN, Niter_CEN]=ChanEsedogluNikolova(Images{i}, double(Masks{i}), lambda_CEN, mu, tho_u, eps, stop_1, stop_2, c1, c2, cinconnu, itermax);
    
    figure(f1);
    subplot(2,NbImages,iter)
    imagesc(Images{i}); axis off; axis image;
    colormap gray
    hold on
    contour(Ub_CP,'r','Linewidth',3);
    hold off
    if iter==ceil(NbImages/2)
        title({'Segmentation par Chambol - Pock';['\lambda= ',num2str(lambda_CP),' , niter= ',num2str(Niter_CP)]});
    else
        title(['\lambda= ',num2str(lambda_CP),' , niter= ',num2str(Niter_CP)]);
    end
    
    subplot(2,NbImages,iter+NbImages)
    imagesc(Images{i}); axis off; axis image;
    colormap gray
    hold on
    contour(Ub_CEN,'r','Linewidth',3);
    hold off
    if iter==ceil(NbImages/2)
        title({'Segmentation par Chan - Esedoglu - Nikolova';['\lambda= ',num2str(lambda_CEN),' , niter= ',num2str(Niter_CEN)]});
    else
        title(['\lambda= ',num2str(lambda_CEN),' , niter= ',num2str(Niter_CEN)]);
    end
    
    figure(f2);
    subplot(3,NbImages,iter);
    plot(1:Niter_CP,J_CP)
    title({['Image ',num2str(i)],['Fonction cout pour \lambda= ',num2str(lambda_CP)]})
    xlabel('iterations')
    subplot(3,NbImages,iter+NbImages)
    semilogy(1:Niter_CP,Err_u_Cp)
    title('Erreur relative sur u_{n+1}-u_n')
    xlabel('iterations')
    subplot(3,NbImages,iter+2*NbImages)
    semilogy(1:Niter_CP,Err_J_CP)
    title('Erreur relative sur J(u_{n+1})-J(u_n)')
    xlabel('iterations')
    
    figure(f3);
    subplot(3,NbImages,iter);
    plot(1:Niter_CEN,J_CEN)
    title({['Image ',num2str(i)],['Fonction cout pour \lambda= ',num2str(lambda_CEN)]})
    xlabel('iterations')
    subplot(3,NbImages,iter+NbImages)
    semilogy(1:Niter_CEN,Err_u_CEN)
    title('Erreur relative sur u_{n+1}-u_n')
    xlabel('iterations')
    subplot(3,NbImages,iter+2*NbImages)
    semilogy(1:Niter_CEN,Err_J_CEN)
    title('Erreur relative sur J(u_{n+1})-J(u_n)')
    xlabel('iterations')
    
% Importation des masques ou creation de ces derniers
if ~Import
    for m=Mask2Change
        Masks{m}=roipoly(Images{m});
    end
    close;
    for b=Background2Change
        Backgrounds{b}=roipoly(Images{b});
    end
    close;
end
%     xlabel('iterations')
end

if ~Import
    save('Images/Masks','Masks');
    save('Images/Backgrounds','Backgrounds');
    save('Images/Regions','Regions');
end