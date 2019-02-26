addpath(genpath('../../'));
Image1=zeros(248);
Image1(62:186,62:186)=255;
Image2=double(imread('eight.tif'));
Image3=double(imread('TumeurCerveaubis.png'));
Image3=Image3(:,:,1);

Images={Image1,Image2,Image3};

%% Ajout de bruit aux images
bruitage=false;
if bruitage==true
    bruit=0.05;
    
    for i=1:3
        Images{i}=Images{i}+bruit*255*randn(size(Images{i}));
    
        Min=min(min(min(Images{i})),0);
        Max=max(max(max(Images{i})),255);

        Images{i}=Images{i}-Min;
        Images{i}=Images{i}*255/(Max-Min);
    end
end

mask1=roipoly(Images{1}/255.);
mask2=roipoly(Images{2}/255.);
mask3=roipoly(Images{3}/255.); close;
Masks={mask1,mask2,mask3};

%% Paramètres numériques
itermax=200; cinconnu=false;
mu=0.1; eps=1;
theta=1;
tho_u=0.5; tho_z=0.25;
stop_1=1.e-6; stop_2=1.e-6;

%% Boucle des lambdas
pow_min=-4; pow_max=2;
      
for ind=pow_min:pow_max
% for ind=1:9
%     lambda=ind*10^-4;
    lambda=10^ind;
    f1=figure('Name','Resultat de la segmentation','NumberTitle','off'); 
    f2=figure('Name','Evolution des quantités de contrôle','NumberTitle','off');
    
    for i=1:3
        c1=mean(Images{i}(Masks{i})); c2=mean(Images{i}(~Masks{i}));
        [Ub, J, Err_u,Err_J, Niter]=DualFormulation(Images{i}, double(Masks{i}), lambda, mu, tho_u, theta, tho_z,stop_1,stop_2, c1, c2, cinconnu, itermax);
        %[Ub, J, Err_u,Err_J, Niter]=ChanEsedogluNikolova(Images{i}, double(Masks{i}), lambda, mu, tho_u, eps, stop_1, stop_2, c1, c2, cinconnu, itermax);

        figure(f1);
        subplot(1,3,i)
        imagesc(Images{i}); axis off; axis image;
        colormap gray
        hold on
        contour(Ub,'r','Linewidth',3);
        title(['Segmentation par C-P avec \lambda= ',num2str(lambda),' et pour ',num2str(Niter),' iterations']);
        hold off
        
        figure(f2);
        subplot(3,3,i);
        plot(1:Niter,J)
        title(['Fonction cout pour \lambda= ',num2str(lambda)])
        xlabel('iterations')
        subplot(3,3,i+3)
        semilogy(1:Niter,Err_u)
        title('Erreur relative sur u_{n+1}-u_n')
        xlabel('iterations')
        subplot(3,3,i+6)
        semilogy(1:Niter,Err_J)
        title('Erreur relative sur J(u_{n+1})-J(u_n)')
        xlabel('iterations')
    end
end