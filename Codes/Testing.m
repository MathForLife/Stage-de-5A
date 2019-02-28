addpath('../')

%% Valeur de lambda à tester + var booléennes
lambda=61;
Import=true;
Bruitage=false;

%% Lecture des images
Image1=zeros(248);
Image1(62:186,62:186)=255;
Image2=double(imread('eight.tif'));
Image3=double(imread('TumeurCerveaubis.png'));

Images={Image1,Image2,Image3};
%% Normalisation des images
for i=1:3
    Images{i}=Image_Normalisation(Images{i},"2D");
end

%% Ajout de bruit aux images
if Bruitage==true
    bruit=0.05;
    
    for i=1:3
        Images{i}=Images{i}+bruit*randn(size(Images{i}));
        Images{i}=Image_Normalisation(Images{i},"2D");
    end
end

% Importation des masques ou création de ces derniers
if Import
    load('Images/Masks','Masks');
else
    mask1=roipoly(Images{1});
    mask2=roipoly(Images{2});
    mask3=roipoly(Images{3}); close;
    Masks={mask1,mask2,mask3};
end

%% Paramètres numériques
itermax=200; cinconnu=false;
mu=0.1; eps=1;
theta=1;
tho_u=0.5; tho_z=0.25;
stop_1=1.e-6; stop_2=1.e-6;

f1=figure('Name','Resultat de la segmentation','NumberTitle','off');
f2=figure('Name','Evolution des quantités de contrôle pour C-P','NumberTitle','off');
f3=figure('Name','Evolution des quantités de contrôle pour C-E-N','NumberTitle','off');
for i=1:3
    c1=mean(Images{i}(Masks{i})); c2=mean(Images{i}(~Masks{i}));
    [Ub_CP, J_CP, Err_u_Cp,Err_J_CP, Niter_CP]=DualFormulation(Images{i}, double(Masks{i}), lambda, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax);
    [Ub_CEN, J_CEN, Err_u_CEN,Err_J_CEN, Niter_CEN]=ChanEsedogluNikolova(Images{i}, double(Masks{i}), lambda, mu, tho_u, eps, stop_1, stop_2, c1, c2, cinconnu, itermax);
    
    figure(f1);
    subplot(2,3,i)
    imagesc(Images{i}); axis off; axis image;
    colormap gray
    hold on
    contour(Ub_CP,'r','Linewidth',3);
    hold off
    if i==2
        title({'Segmentation par Chambol-Pock';['\lambda= ',num2str(lambda),' , niter= ',num2str(Niter_CP)]});
    else
        title(['\lambda= ',num2str(lambda),' , niter= ',num2str(Niter_CP)]);
    end
    
    subplot(2,3,i+3)
    imagesc(Images{i}); axis off; axis image;
    colormap gray
    hold on
    contour(Ub_CEN,'r','Linewidth',3);
    hold off
    if i==2
        title({'Segmentation par Chan-Esedoglu-Nikolova';['\lambda= ',num2str(lambda),' , niter= ',num2str(Niter_CEN)]});
    else
        title(['\lambda= ',num2str(lambda),' , niter= ',num2str(Niter_CEN)]);
    end
    
    figure(f2);
    subplot(3,3,i);
    plot(1:Niter_CP,J_CP)
    title(['Fonction cout pour \lambda= ',num2str(lambda)])
    xlabel('iterations')
    subplot(3,3,i+3)
    semilogy(1:Niter_CP,Err_u_Cp)
    title('Erreur relative sur u_{n+1}-u_n')
    xlabel('iterations')
    subplot(3,3,i+6)
    semilogy(1:Niter_CP,Err_J_CP)
    title('Erreur relative sur J(u_{n+1})-J(u_n)')
    xlabel('iterations')
    
    figure(f3);
    subplot(3,3,i);
    plot(1:Niter_CEN,J_CEN)
    title(['Fonction cout pour \lambda= ',num2str(lambda)])
    xlabel('iterations')
    subplot(3,3,i+3)
    semilogy(1:Niter_CEN,Err_u_CEN)
    title('Erreur relative sur u_{n+1}-u_n')
    xlabel('iterations')
    subplot(3,3,i+6)
    semilogy(1:Niter_CEN,Err_J_CEN)
    title('Erreur relative sur J(u_{n+1})-J(u_n)')
    xlabel('iterations')
end

if ~Import
    save('Images/Masks','Masks');
end